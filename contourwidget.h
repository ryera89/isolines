#ifndef CONTOURWIDGET_H
#define CONTOURWIDGET_H

#include <QWidget>
#include <QPair>
#include "marching_squares.h"

class PlotSettings
{
public:
    PlotSettings();
    void scroll(int dx, int dy);
    void adjust();
    double spanX() const { return maxX - minX; }
    double spanY() const { return maxY - minY; }
    double minX;
    double maxX;
    int numXTicks;
    double minY;
    double maxY;
    int numYTicks;
private:
    static void adjustAxis(double &min, double &max, int &numTicks);
};


class contourWidget : public QWidget
{
    Q_OBJECT
public:
    explicit contourWidget(QWidget *parent = nullptr);
    void contour(const VecDoub &x,const VecDoub &y, const MatDoub &Z,const double &isovalue);
    void contour(const VecDoub &x,const VecDoub &y, const MatDoub &Z,const QVector<double> &isovalues);

signals:

protected:
    void paintEvent(QPaintEvent *event);
    QSize sizeHint() const;

    void drawGrid(QPainter *painter);
    void drawIsoline(QPainter *painter);

public slots:
    void addIsoline(const QVector<QPair<QPointF,QPointF>> &isoline);
    void addIsolines(const QVector<QVector<QPair<QPointF,QPointF>>> &isolines);
    void clearIsolines();
private:
    PlotSettings settings;
    const int MARGIN = 50;
    QVector<QVector<QPair<QPointF,QPointF>>> v_isolines;
    marching_squares mSquares;

};

#endif // CONTOURWIDGET_H
