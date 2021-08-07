#ifndef ISOLINES_WIDGET_H
#define ISOLINES_WIDGET_H

#include <QWidget>
#include <QtCharts/QtCharts>
#include "linalg/mathmatrix.h"

using namespace QtCharts;

class isolines_widget : public QWidget
{
    Q_OBJECT
public:
    explicit isolines_widget(QWidget *parent = nullptr);

    void drawContours(const MatDoub &X,const MatDoub &Y, const MatDoub &Z);
    void addContour(const QList<QPointF> &contourPath);
    void addConstraint(const QList<QPointF> &constPath);
    void addConvergenceLine(const QList<QPointF> &convPath);
    QChartView* chartView(){return view;}

signals:

public slots:
private:
    QChartView *view;
    QChart *chart;
    QVector<QLineSeries*> v_countour; //contornos
    QVector<QLineSeries*> v_constrains; //restricciones
    QVector<QLineSeries*> v_csteeps; //convergencia

};

#endif // ISOLINES_WIDGET_H
