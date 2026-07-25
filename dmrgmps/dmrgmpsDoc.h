
// dmrgDoc.h : interface of the CdmrgmpsDoc class
//


#pragma once

#include "Chart.h"

#include "ComputationThread.h"
#include "Options.h"

class CdmrgmpsView;

class CdmrgmpsDoc : public CDocument
{
protected: // create from serialization only
	CdmrgmpsDoc();
	DECLARE_DYNCREATE(CdmrgmpsDoc)

// Attributes
public:
	Chart m_Chart;
// Operations
	bool IsFinished();
	void UpdateChartData();
	CdmrgmpsView* GetView();
	void StartComputing();
	void SetYAxisRange();
// Overrides
	BOOL OnNewDocument() override;
	void Serialize(CArchive& ar) override;
#ifdef SHARED_HANDLERS
	void InitializeSearchContent() override;
	void OnDrawThumbnail(CDC& dc, LPRECT lprcBounds) override;
#endif // SHARED_HANDLERS

// Implementation
	~CdmrgmpsDoc() override;
#ifdef _DEBUG
	void AssertValid() const override;
	void Dump(CDumpContext& dc) const override;
#endif

private:
	ComputationThread *thread;
	Options options;

// Generated message map functions
	DECLARE_MESSAGE_MAP()

#ifdef SHARED_HANDLERS
	// Helper function that sets search content for a Search Handler
	void SetSearchContent(const CString& value);
#endif // SHARED_HANDLERS
};
