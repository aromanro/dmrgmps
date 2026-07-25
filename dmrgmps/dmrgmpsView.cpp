
// dmrgView.cpp : implementation of the CdmrgmpsView class
//

#include "stdafx.h"
// SHARED_HANDLERS can be defined in an ATL project implementing preview, thumbnail
// and search filter handlers and allows sharing of document code with that project.
#ifndef SHARED_HANDLERS
#include "dmrgmps.h"
#endif

#include "dmrgmpsDoc.h"
#include "dmrgmpsView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CdmrgmpsView

IMPLEMENT_DYNCREATE(CdmrgmpsView, CView)

BEGIN_MESSAGE_MAP(CdmrgmpsView, CView)
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CdmrgmpsView::OnFilePrintPreview)
	ON_WM_CONTEXTMENU()
	ON_WM_RBUTTONUP()
	ON_WM_ERASEBKGND()
	ON_WM_SETCURSOR()
	ON_WM_DESTROY()
	ON_WM_TIMER()
END_MESSAGE_MAP()


BOOL CdmrgmpsView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

// CdmrgmpsView drawing

void CdmrgmpsView::OnDraw(CDC* pDC)
{
	CdmrgmpsDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO: add draw code for native data here

	CRect rect;
	GetClientRect(rect);

	pDoc->m_Chart.Draw(pDC, rect);
}


// CdmrgmpsView printing


void CdmrgmpsView::OnFilePrintPreview()
{
#ifndef SHARED_HANDLERS
	AFXPrintPreview(this);
#endif
}

BOOL CdmrgmpsView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CdmrgmpsView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add extra initialization before printing
}

void CdmrgmpsView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}

void CdmrgmpsView::OnRButtonUp(UINT /* nFlags */, CPoint point)
{
	ClientToScreen(&point);
	OnContextMenu(this, point);
}

void CdmrgmpsView::OnContextMenu(CWnd* /* pWnd */, CPoint /*point*/)
{
#ifndef SHARED_HANDLERS
#endif
}


// CdmrgmpsView diagnostics

#ifdef _DEBUG
void CdmrgmpsView::AssertValid() const
{
	CView::AssertValid();
}

void CdmrgmpsView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CdmrgmpsDoc* CdmrgmpsView::GetDocument() const // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CdmrgmpsDoc)));
	return dynamic_cast<CdmrgmpsDoc*>(m_pDocument);
}
#endif //_DEBUG


// CdmrgmpsView message handlers


void CdmrgmpsView::OnPrepareDC(CDC* pDC, CPrintInfo* pInfo)
{
	CView::OnPrepareDC(pDC, pInfo);

	if (pDC->IsPrinting())
	{
		CRect rect;
		GetClientRect(rect);

		pDC->SetMapMode(MM_ISOTROPIC);

		int cx = pDC->GetDeviceCaps(HORZRES);
		int cy = pDC->GetDeviceCaps(VERTRES);
		
		pDC->SetWindowExt(rect.Width(), rect.Height());
		pDC->SetViewportExt(cx, cy);
		pDC->SetViewportOrg(0, 0);
	}
}


BOOL CdmrgmpsView::OnEraseBkgnd(CDC* /*pDC*/)
{
	return TRUE;
}


BOOL CdmrgmpsView::OnSetCursor(CWnd* pWnd, UINT nHitTest, UINT message)
{
	CdmrgmpsDoc* pDoc = GetDocument();

	if (pDoc && !pDoc->IsFinished())
	{
		RestoreWaitCursor();

		return TRUE;
	}

	return CView::OnSetCursor(pWnd, nHitTest, message);
}


void CdmrgmpsView::OnDestroy()
{
	CView::OnDestroy();

	if (timer) KillTimer(timer);
}


void CdmrgmpsView::OnTimer(UINT_PTR nIDEvent)
{
	CView::OnTimer(nIDEvent);

	CdmrgmpsDoc* pDoc = GetDocument();

	if (pDoc->IsFinished())
	{
		KillTimer(timer);
		timer = 0;

		pDoc->UpdateChartData();

		EndWaitCursor();
		Invalidate();
	}	
}


void CdmrgmpsView::StartTimer()
{
	if (!timer) timer = SetTimer(1, 1000, nullptr);
	BeginWaitCursor();
}
