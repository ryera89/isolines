#include "isolines_widget.h"
#include <QHBoxLayout>
#include <QHash>
isolines_widget::isolines_widget(QWidget *parent) : QWidget(parent)
{   
    chart = new QChart;
    view = new QChartView(chart);
    view->chart()->setTheme(QChart::ChartThemeDark);
    view->setRenderHint(QPainter::Antialiasing);
    view->chart()->legend()->setVisible(false);
    view->chart()->createDefaultAxes();

    QHBoxLayout *layout = new QHBoxLayout;
    layout->addWidget(view);
    this->setLayout(layout);
}

void isolines_widget::drawContours(const MatDoub &X, const MatDoub &Y, const MatDoub &Z)
{
    assert(X.rows() == Y.rows() && X.cols() == Y.cols() && X.rows() == Z.rows() && X.cols() == Z.cols());

    QHash<double,QList<QPointF>> grid;

    for (size_t i = 0; i < Z.rows(); ++i){
        for (size_t j = 0; j < Z.cols(); ++j){
            grid[Z(i,j)].append(QPointF(X(i,j),Y(i,j)));
        }
    }
    int count = 0;
    for (QHash<double,QList<QPointF>>::iterator its = grid.begin(); its != grid.end(); ++its){
        if (its.value().size() <= 6) ++count;
    }
    for (QHash<double,QList<QPointF>>::iterator its = grid.begin(); its != grid.end(); ++its){
        if (its.value().size() > 6){
            addContour(its.value());
        }
    }
}
void isolines_widget::addContour(const QList<QPointF> &contourPath)
{
    QLineSeries *serie = new QLineSeries;
    serie->append(contourPath);
    serie->setColor(Qt::green);

    //QLineSeries *serie2 = new QSplineSeries;
    //serie2->append(isoline2);
    //QAreaSeries *serie = new QAreaSeries(serie1,serie2);
    //serie->setOpacity(0.2);
    //QColor color = serie->color();
    //QBrush brush(Qt::transparent);
    //serie->setColor(Qt::transparent);
    //QPen pen = serie->pen();
    //pen.setWidth(3);
    //serie->setPen(pen);
    //serie->setBorderColor(Qt::darkGreen);
    v_countour.push_back(serie);
    view->chart()->addSeries(serie);
    view->chart()->createDefaultAxes();
    //view->chart()->axisX()->setRange(-5,5);
    //view->chart()->axisY()->setRange(-5,5);
}

void isolines_widget::addConstraint(const QList<QPointF> &constPath)
{
    QLineSeries *serie = new QLineSeries;
    serie->append(constPath);
    serie->setColor(Qt::darkRed);

    v_constrains.push_back(serie);
    view->chart()->addSeries(serie);
    view->chart()->createDefaultAxes();
}

void isolines_widget::addConvergenceLine(const QList<QPointF> &convPath)
{
    QLineSeries *serie = new QLineSeries;
    serie->append(convPath);

    v_csteeps.push_back(serie);
}
