#include "MeshParamWidget.h"

MeshParamWidget::MeshParamWidget(QWidget *parent)
	: QWidget(parent)
{
	CreateTabWidget();
	CreateLayout();
}

MeshParamWidget::~MeshParamWidget()
{
}

void MeshParamWidget::CreateTabWidget(void)
{
	pbPrintInfo = new QPushButton(tr("Print Mesh Info"));
	connect(pbPrintInfo, SIGNAL(clicked()), SIGNAL(PrintInfoSignal()));

	pMarkUVDraw = new QCheckBox(tr("Draw Mesh UV"));
	connect(pMarkUVDraw, SIGNAL(clicked()), SIGNAL(ClickDrawUVSignal()));

	pUVScaleBox = new QSpinBox();
	pUVScaleBox->setMinimum(0);
	pUVScaleBox->setMaximum(100);
	pUVScaleBox->setSingleStep(10);

	pUVScale = new QSlider(Qt::Horizontal);
	pUVScale->setMinimum(0);
	pUVScale->setMaximum(100);
	pUVScale->setSingleStep(10);
	pUVScale->setTickPosition(QSlider::TicksAbove);
	connect(pUVScaleBox, SIGNAL(valueChanged(int)), pUVScale, SLOT(setValue(int)));
	connect(pUVScale, SIGNAL(valueChanged(int)), pUVScaleBox, SLOT(setValue(int)));
	
	connect(pUVScale, SIGNAL(valueChanged(int)), SIGNAL(ChangeUVScaleSignal(int)));

	QHBoxLayout* hLayout = new QHBoxLayout();
	hLayout->addWidget(pUVScaleBox);
	hLayout->addWidget(pUVScale);

	QVBoxLayout *layout = new QVBoxLayout();
	layout->addWidget(pbPrintInfo);
	layout->addWidget(pMarkUVDraw);
	//layout->addWidget(pUVScale);
	//layout->addWidget(hLayout);
	layout->addLayout(hLayout);
	layout->addStretch();
	wParam = new QWidget();
	wParam->setLayout(layout);
	saParam = new QScrollArea();
	saParam->setFocusPolicy(Qt::NoFocus);
	saParam->setFrameStyle(QFrame::NoFrame);
	saParam->setWidget(wParam);
	saParam->setWidgetResizable(true);
}

void MeshParamWidget::CreateLayout(void)
{
	twParam = new QTabWidget();
	twParam->addTab(saParam, "Tab");
	QGridLayout *layout = new QGridLayout();
	layout->addWidget(twParam, 0, 0, 1, 1);
	this->setLayout(layout);
}