#include "isolineswindow.h"
#include <QPushButton>
#include <QComboBox>
#include <QTextEdit>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QGroupBox>
#include <QGridLayout>
#include <QLabel>
#include <QDoubleSpinBox>
#include "isolineslineedit.h"
#include <QMenu>
#include <QAction>

//Funciones a graficar isolineas
#include "examples/ejem_16_04_interior_point.h"
#include "examples/ejemplo_rosenbrock_interior_point.h"

template<typename T>
bool equal(const T &v1, const T &v2){
    return v1 == v2;
}
bool applyDiff(const QVector<double> &v,const double &val,bool (*pred)(const double &,const double &)){
    bool flag = true;
    for (int i = 0; i < v.size(); ++i){
        if (pred(v[i],val)) return false;
    }
    return flag;
}


isolinesWindow::isolinesWindow(QWidget *parent) : QMainWindow(parent)
{
    m_contourWidget = new contourWidget;
    m_contourWidget->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);

    m_resultsDisplay = new QTextEdit;
    m_resultsDisplay->setSizePolicy(QSizePolicy::Minimum,QSizePolicy::Expanding);

    m_optTypeComboBox = new QComboBox;
    m_optTypeComboBox->addItem("Optimizacion sin restricciones");
    m_optTypeComboBox->addItem("Optimizacion con restricciones");

    m_withOutConstGroupBox = new QGroupBox(tr("Sin Restricciones"));
    m_problemWORComboBox = new QComboBox;
    m_problemWORComboBox->addItem(tr("Funcion de Rosenbrock"));
    m_methodWORComboBox = new QComboBox;
    m_methodWORComboBox->addItem(tr("Steepest Descent"));
    m_methodWORComboBox->addItem(tr("Gradiente Conjugado FR-PR"));
    m_methodWORComboBox->addItem(tr("Gradiente Conjugado Modificado"));
    m_methodWORComboBox->addItem(tr("Quasi-Newton BFGS"));


    QVBoxLayout *worlayout = new QVBoxLayout;
    worlayout->addWidget(m_problemWORComboBox);
    worlayout->addWidget(m_methodWORComboBox);

    m_withOutConstGroupBox->setLayout(worlayout);

    m_withConstGroupBox = new QGroupBox(tr("Con Restricciones"));
    m_problemWRComboBox = new QComboBox;
    m_methodWRComboBox = new QComboBox;

    //m_qp_eqConstAction = new QAction(tr("QP Restricciones de igualdad"));
    //m_qp_actSetAction = new QAction(tr("QP Set Activo"));
    //m_qp_gradProjAction = new QAction(tr("QP Gradiente Proyectado"));
    //m_interiorPointAction = new QAction(tr("Punto Interior"));

    //m_methodWRComboBox->addAction(m_qp_eqConstAction);
    //m_methodWRComboBox->addAction(m_qp_actSetAction);
    //m_methodWRComboBox->addAction(m_qp_gradProjAction);
    //m_methodWRComboBox->addAction(m_interiorPointAction);


    QVBoxLayout *wrlayout = new QVBoxLayout;
    wrlayout->addWidget(m_problemWRComboBox);
    wrlayout->addWidget(m_methodWRComboBox);

    m_withConstGroupBox->setLayout(wrlayout);

    m_minimizePushButton = new QPushButton(tr("Minimizar"));

    QWidget *centralwidget = new QWidget;

    QHBoxLayout *clayout = new QHBoxLayout;
    clayout->addWidget(m_contourWidget);
    clayout->addWidget(m_resultsDisplay);

    m_isolinesGroupBox = new QGroupBox(tr("Isolineas"));
    m_XLabel = new QLabel(tr("X"));
    m_YLabel = new QLabel(tr("Y"));
    m_rangoLabel = new QLabel(tr("Rango"));
    m_stepLabel = new QLabel(tr("Paso"));
    m_isolinesLabel = new QLabel(tr("Z"));

    m_xbSpinButton = new QDoubleSpinBox;
    m_xbSpinButton->setRange(-10000,10000);
    m_xbSpinButton->setValue(0.0);
    m_xeSpinButton = new QDoubleSpinBox;
    m_xeSpinButton->setValue(10.0);
    m_xeSpinButton->setRange(-10000,10000);
    m_ybSpinButton = new QDoubleSpinBox;
    m_ybSpinButton->setValue(0.0);
    m_ybSpinButton->setRange(-10000,10000);
    m_yeSpinButton = new QDoubleSpinBox;
    m_yeSpinButton->setValue(10.0);
    m_yeSpinButton->setRange(-10000,10000);
    m_stepSpinBox = new QDoubleSpinBox;
    m_stepSpinBox->setValue(0.1);
    m_stepSpinBox->setMinimum(0.05);
    m_stepSpinBox->setSingleStep(0.05);

    m_isolinesLineEdit = new isolinesLineEdit;

    m_obtenerIsolineasAction = new QAction(tr("Graficar nuevas isolineas"));
    connect(m_obtenerIsolineasAction,SIGNAL(triggered(bool)),this,SLOT(graficarNuevasIsolineas()));
    m_agregarIsolineasAction = new QAction(tr("Agregar nuevas isolineas"));
    connect(m_agregarIsolineasAction,SIGNAL(triggered(bool)),this,SLOT(agregarNuevasIsolineas()));
    m_isolineasMenu = new QMenu;
    m_isolineasMenu->addAction(m_obtenerIsolineasAction);
    m_isolineasMenu->addAction(m_agregarIsolineasAction);

    m_getIsolinesButton = new QPushButton(tr("Graficar Isolineas"));
    m_getIsolinesButton->setMenu(m_isolineasMenu);

    QGridLayout *rangeLayout = new QGridLayout;
    rangeLayout->addWidget(m_rangoLabel,0,1,1,2,Qt::AlignCenter);
    rangeLayout->addWidget(m_XLabel,1,0,1,1);
    rangeLayout->addWidget(m_YLabel,2,0,1,1);
    rangeLayout->addWidget(m_xbSpinButton,1,1,1,1);
    rangeLayout->addWidget(m_xeSpinButton,1,2,1,1);
    rangeLayout->addWidget(m_ybSpinButton,2,1,1,1);
    rangeLayout->addWidget(m_yeSpinButton,2,2,1,1);
    rangeLayout->addWidget(m_stepLabel,3,0,1,1);
    rangeLayout->addWidget(m_stepSpinBox,3,1,1,2);
    rangeLayout->addWidget(m_isolinesLabel,4,0,1,1);
    rangeLayout->addWidget(m_isolinesLineEdit,4,1,1,2);
    rangeLayout->addWidget(m_getIsolinesButton,5,0,1,3);
    m_isolinesGroupBox->setLayout(rangeLayout);

    QVBoxLayout *rlayout = new QVBoxLayout;
    rlayout->addWidget(m_optTypeComboBox);
    rlayout->addWidget(m_withOutConstGroupBox);
    rlayout->addWidget(m_withConstGroupBox);
    rlayout->addWidget(m_isolinesGroupBox);
    rlayout->addWidget(m_minimizePushButton);
    rlayout->addStretch();

    clayout->addLayout(rlayout);

    centralwidget->setLayout(clayout);
    setCentralWidget(centralwidget);
}

void isolinesWindow::graficarNuevasIsolineas()
{
    //TODO aqui va agragar por la eleccion que de haya hecho
    int optType = m_optTypeComboBox->currentIndex();
    if (optType == 0){
        ejemplo_rosenbrock func;
        m_contourWidget->clearIsolines();
        m_isovalues.clear();
        this->computeIsolines(func);
    }else{
        ejemplo_16_04 func;
        m_contourWidget->clearIsolines();
        m_isovalues.clear();
        this->computeIsolines(func);
    }

}

void isolinesWindow::agregarNuevasIsolineas()
{
    int optType = m_optTypeComboBox->currentIndex();
    if (optType == 0){
        ejemplo_rosenbrock func;
        this->computeIsolines(func);
    }else{
        ejemplo_16_04 func;
        this->computeIsolines(func);
    }
}

template<typename T>
void isolinesWindow::computeIsolines(T &func)
{
    double xb = m_xbSpinButton->value();
    double xe = m_xeSpinButton->value();
    double yb = m_ybSpinButton->value();
    double ye = m_yeSpinButton->value();

    double step = m_stepSpinBox->value();

    double xrange = xe - xb;
    double yrange = ye - yb;

    int xsize = std::floor(xrange/step) + 1;
    int ysize = std::floor(yrange/step) + 1;

    VecDoub x(xsize);
    VecDoub y(ysize);

    for (int i = 0; i < xsize; ++i) x[i] = xb + step*i;
    for (int i = 0; i < ysize; ++i) y[i] = yb + step*i;

    MatDoub Z(ysize,xsize);

    VecDoub xi(2);

    for (int i = 0; i < ysize; ++i){
        for (int j = 0; j < xsize; ++j){
            xi[0] = x[j];
            xi[1] = y[i];
            Z(i,j) = func(xi);
        }
    }
    QVector<double> isovalues = m_isolinesLineEdit->getIsolinesValues();

    if (!m_isovalues.size()){
        m_isovalues.append(isovalues);
    }else{
        for (int i = 0; i < isovalues.size(); ++i){
            if (applyDiff(m_isovalues,isovalues[i],equal)) m_isovalues.append(isovalues[i]);
        }
    }

    m_contourWidget->contour(x,y,Z,m_isovalues);
}
