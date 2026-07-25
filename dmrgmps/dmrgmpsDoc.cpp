
// dmrgDoc.cpp : implementation of the CdmrgmpsDoc class
//

#include "stdafx.h"
// SHARED_HANDLERS can be defined in an ATL project implementing preview, thumbnail
// and search filter handlers and allows sharing of document code with that project.
#ifndef SHARED_HANDLERS
#include "dmrgmps.h"
#endif

#include "dmrgmpsDoc.h"
#include "dmrgmpsView.h"

#include "DMRGHeisenbergSpinOneHalf.h"
#include "DMRGHeisenbergSpinOne.h"

#include "DMRGThread.h"

#include <propkey.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CdmrgmpsDoc

IMPLEMENT_DYNCREATE(CdmrgmpsDoc, CDocument)

BEGIN_MESSAGE_MAP(CdmrgmpsDoc, CDocument)
END_MESSAGE_MAP()


// CdmrgmpsDoc construction/destruction

CdmrgmpsDoc::CdmrgmpsDoc()
	: thread(nullptr)
{
	// TODO: add one-time construction code here
	m_Chart.useSpline = false;
	m_Chart.XAxisLabel = L"Site";
	m_Chart.YAxisLabel = L"<S_i S_{i+1}>";

	m_Chart.XAxisMin = 0;
}

CdmrgmpsDoc::~CdmrgmpsDoc()
{
	delete thread;
}

BOOL CdmrgmpsDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	SetTitle(L"No data");

	return TRUE;
}




// CdmrgmpsDoc serialization

void CdmrgmpsDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}

#ifdef SHARED_HANDLERS

// Support for thumbnails
void CdmrgmpsDoc::OnDrawThumbnail(CDC& dc, LPRECT lprcBounds)
{
	// Modify this code to draw the document's data
	dc.FillSolidRect(lprcBounds, RGB(255, 255, 255));

	CString strText = _T("TODO: implement thumbnail drawing here");
	LOGFONT lf;

	CFont* pDefaultGUIFont = CFont::FromHandle((HFONT) GetStockObject(DEFAULT_GUI_FONT));
	pDefaultGUIFont->GetLogFont(&lf);
	lf.lfHeight = 36;

	CFont fontDraw;
	fontDraw.CreateFontIndirect(&lf);

	CFont* pOldFont = dc.SelectObject(&fontDraw);
	dc.DrawText(strText, lprcBounds, DT_CENTER | DT_WORDBREAK);
	dc.SelectObject(pOldFont);
}

// Support for Search Handlers
void CdmrgmpsDoc::InitializeSearchContent()
{
	CString strSearchContent;
	// Set search contents from document's data. 
	// The content parts should be separated by ";"

	// For example:  strSearchContent = _T("point;rectangle;circle;ole object;");
	SetSearchContent(strSearchContent);
}

void CdmrgmpsDoc::SetSearchContent(const CString& value)
{
	if (value.IsEmpty())
	{
		RemoveChunk(PKEY_Search_Contents.fmtid, PKEY_Search_Contents.pid);
	}
	else
	{
		CMFCFilterChunkValueImpl *pChunk = nullptr;
		ATLTRY(pChunk = new CMFCFilterChunkValueImpl);
		if (pChunk != nullptr)
		{
			pChunk->SetTextValue(PKEY_Search_Contents, value, CHUNK_TEXT);
			SetChunkValue(pChunk);
		}
	}
}

#endif // SHARED_HANDLERS

// CdmrgmpsDoc diagnostics

#ifdef _DEBUG
void CdmrgmpsDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CdmrgmpsDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG


// CdmrgmpsDoc commands


bool CdmrgmpsDoc::IsFinished()
{
	return !thread || thread->terminated;
}


void CdmrgmpsDoc::UpdateChartData()
{
	if (!thread || !thread->terminated) return;

	std::list<double> results;

	CString spinStr;
	if (0 == options.model)
	{
		results.swap(((DMRGThread<DMRG::Heisenberg::DMRGHeisenbergSpinOneHalf>*)thread)->dmrg.results);
		spinStr = L"1/2";
	}
	else
	{
		results.swap(((DMRGThread<DMRG::Heisenberg::DMRGHeisenbergSpinOne>*)thread)->dmrg.results);
		spinStr = L"1";
	}

	double result = thread->result;
	double gapResult = thread->gapResult;

	delete thread;
	thread = nullptr;

	int chainLength = options.sites;
	SetTitle(L"Finished");

	std::vector<double> x, y;
	int i = 0;
	for (auto res : results)
	{
		x.push_back(++i);
		y.push_back(res);
	}

	m_Chart.clear();

	m_Chart.AddDataSet(x.data(), y.data(), static_cast<int>(x.size()), 2, RGB(255,0,0));

	m_Chart.SetNumBigTicksX(theApp.options.bigTicksX);
	m_Chart.SetNumBigTicksY(theApp.options.bigTicksY);
	m_Chart.SetNumTicksX(theApp.options.bigTicksX * theApp.options.smallTicksX);
	m_Chart.SetNumTicksY(theApp.options.bigTicksY * theApp.options.smallTicksY);


	if (options.calculateEnergyGap)
		m_Chart.title.Format(L"L=%d, S=%s, Open BCs, Energy/site: %.5f, Energy Gap: %.5f", chainLength, (const wchar_t*)spinStr, result / chainLength, gapResult);	
	else
		m_Chart.title.Format(L"L=%d, S=%s, Open BCs, Energy/site: %.5f", chainLength, (const wchar_t*)spinStr, result / chainLength);	

	m_Chart.XAxisMax = static_cast<int>(x.size() + 1);

	SetYAxisRange();
}


CdmrgmpsView* CdmrgmpsDoc::GetView()
{
	POSITION pos = GetFirstViewPosition();
	while (pos)
	{
		CView* pView = GetNextView(pos);
		if (pView->IsKindOf(RUNTIME_CLASS(CdmrgmpsView)))
			return dynamic_cast<CdmrgmpsView*>(pView);
	}

	return nullptr;
}

void CdmrgmpsDoc::StartComputing()
{
	if (thread) return;

	options = theApp.options;

	SetTitle(L"Computing");

	CdmrgmpsView* view = GetView();
	if (view)
	{
		view->StartTimer();			
	}

	if (0 == theApp.options.model)
		thread = new DMRGThread<DMRG::Heisenberg::DMRGHeisenbergSpinOneHalf>(options.sites, options.Jz, options.Jxy, options.sweeps, options.states, options.calculateEnergyGap ? options.nrStates : 0, options.method);
	else
		thread = new DMRGThread<DMRG::Heisenberg::DMRGHeisenbergSpinOne>(options.sites, options.Jz, options.Jxy, options.sweeps, options.states, options.calculateEnergyGap ? options.nrStates : 0, options.method);

	thread->Start();
}


void CdmrgmpsDoc::SetYAxisRange()
{
	m_Chart.YAxisMin = theApp.options.minEnergy;
	m_Chart.YAxisMax = theApp.options.maxEnergy;
}
