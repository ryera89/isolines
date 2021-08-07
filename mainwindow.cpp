#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "QGridLayout"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    rosenbrockIso = new isolines_widget;
    ejemplo_1604Iso = new isolines_widget;
    ejemplo_1611Iso = new isolines_widget;
    ejemplo_juanmaIso = new isolines_widget;

    QWidget *centralWidget = new QWidget;
    //QLabel *label = new QLabel("Hello Qt");
    QGridLayout *layout = new QGridLayout;
    layout->addWidget(rosenbrockIso,0,0);
    layout->addWidget(ejemplo_1604Iso,0,1);
    layout->addWidget(ejemplo_1611Iso,1,0);
    layout->addWidget(ejemplo_juanmaIso,1,1);
    centralWidget->setLayout(layout);
    this->setCentralWidget(centralWidget);



}

MainWindow::~MainWindow()
{
    delete ui;
}



