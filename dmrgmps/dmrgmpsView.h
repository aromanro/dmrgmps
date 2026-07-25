
// dmrgView.h : interface of the CdmrgmpsView class
//

#pragma once


class CdmrgmpsView : public CView
{
protected: // create from serialization only
	CdmrgmpsView() = default;
	DECLARE_DYNCREATE(CdmrgmpsView)

// Attributes
	UINT_PTR timer = 0;
public:
	CdmrgmpsDoc* GetDocument() const;

// Operations
	void StartTimer();


// Implementation
#ifdef _DEBUG
	void AssertValid() const override;
	void Dump(CDumpContext& dc) const override;
#endif

private:
	// Overrides
	void OnDraw(CDC* pDC) override;  // overridden to draw this view
	BOOL PreCreateWindow(CREATESTRUCT& cs) override;
	BOOL OnPreparePrinting(CPrintInfo* pInfo) override;
	void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo) override;
	void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo) override;

// Generated message map functions
	afx_msg void OnFilePrintPreview();
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnContextMenu(CWnd* pWnd, CPoint point);
	DECLARE_MESSAGE_MAP()

	void OnPrepareDC(CDC* pDC, CPrintInfo* pInfo = nullptr) override;
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	afx_msg BOOL OnSetCursor(CWnd* pWnd, UINT nHitTest, UINT message);
	afx_msg void OnDestroy();
	afx_msg void OnTimer(UINT_PTR nIDEvent);
};

#ifndef _DEBUG  // debug version in dmrgView.cpp
inline CdmrgmpsDoc* CdmrgmpsView::GetDocument() const
   { return reinterpret_cast<CdmrgmpsDoc*>(m_pDocument); }
#endif

