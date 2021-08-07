#include "contourwidget.h"
#include <QPen>
#include <QPainter>
#include <cmath>
#include <iostream>
contourWidget::contourWidget(QWidget *parent) : QWidget(parent)
{
    setBackgroundRole(QPalette::Dark);
    setAutoFillBackground(true);

}

void contourWidget::contour(const VecDoub &x, const VecDoub &y, const MatDoub &Z, const double &isovalue)
{
    if (!x.size() || !y.size()) return; //TODO que esto lance un error
    size_t n = x.size();
    size_t m = y.size();
    if (x[0] > x[n-1]){
        settings.minX = x[n-1] ;
        settings.maxX = x[0];
    }else{
        settings.minX = x[0] ;
        settings.maxX = x[n-1];
    }
    if (y[0] > y[m-1]){
        settings.minY = x[m-1] ;
        settings.maxY = x[0];
    }else{
        settings.minY = y[0] ;
        settings.maxY = y[m-1];
    }
    settings.adjust();
    this->addIsoline(mSquares.getContour(x,y,Z,isovalue));
}

void contourWidget::contour(const VecDoub &x, const VecDoub &y, const MatDoub &Z, const QVector<double> &isovalues)
{
    if (!x.size() || !y.size()) return; //TODO que esto lance un error
    size_t n = x.size();
    size_t m = y.size();
    if (x[0] > x[n-1]){
        settings.minX = x[n-1] ;
        settings.maxX = x[0];
    }else{
        settings.minX = x[0] ;
        settings.maxX = x[n-1];
    }
    if (y[0] > y[m-1]){
        settings.minY = x[m-1] ;
        settings.maxY = x[0];
    }else{
        settings.minY = y[0] ;
        settings.maxY = y[m-1];
    }
    settings.adjust();
    QVector<QVector<QPair<QPointF,QPointF>>> isolines;
    isolines.reserve(isovalues.size());
    for (int i = 0; i < isovalues.size(); ++i){
        isolines.append(mSquares.getContour(x,y,Z,isovalues[i]));
    }
    this->addIsolines(isolines);
}

void contourWidget::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    painter.save();
    QPen pen(palette().light().color());
    pen.setWidth(1);
    painter.setPen(pen);
    QRect r = rect().adjusted(MARGIN,MARGIN,-MARGIN,-MARGIN);
    painter.drawRect(r);
    painter.restore();
    painter.save();
    drawGrid(&painter);
    painter.restore();
    painter.save();
    drawIsoline(&painter);
    painter.restore();
}

QSize contourWidget::sizeHint() const
{
    return QSize(700,700);
}

void contourWidget::drawGrid(QPainter *painter)
{
    QRect rect(MARGIN, MARGIN, width() - 2*MARGIN + 1, height() - 2*MARGIN + 1);

    //QPen quiteDark = palette().dark().color().light();
    QPen light = palette().light().color();

    painter->setPen(light);

    for (int i = 0; i <= settings.numXTicks; ++i){
        int x = rect.left() + (i*(rect.width()-1)/settings.numXTicks);
        double label = settings.minX + (i*settings.spanX()/settings.numXTicks);
        label = std::floor(label*100)/100;

        painter->drawLine(x,rect.bottom(),x,rect.bottom()+10);
        painter->drawText(x-50,rect.bottom()+10,100,20,Qt::AlignHCenter | Qt::AlignTop,QString::number(label,'f',1));

        painter->drawLine(x,rect.top()-10,x,rect.top());
        painter->drawText(x-50,rect.top()-30,100,20,Qt::AlignHCenter | Qt::AlignBottom,QString::number(label,'f',1));
    }
    for (int i = 0; i <= settings.numYTicks; ++i){
        int y = rect.bottom() - (i*(rect.height()-1)/settings.numYTicks);

        double label = settings.minY + (i*settings.spanY()/settings.numYTicks);
        label = std::floor(label*100)/100;

        painter->drawLine(rect.left()-10,y,rect.left(),y);
        painter->drawText(rect.left()-MARGIN,y-10,MARGIN-12,20,Qt::AlignVCenter | Qt::AlignRight,QString::number(label,'f',1));

        painter->drawLine(rect.right(),y,rect.right()+10,y);
        painter->drawText(rect.right()+12,y-10,MARGIN,20,Qt::AlignVCenter | Qt::AlignLeft,QString::number(label,'f',1));
    }
}

void contourWidget::drawIsoline(QPainter *painter)
{
    QRect rect(MARGIN,MARGIN,width() - 2*MARGIN,height() - 2*MARGIN);

    painter->setClipRect(rect.adjusted(1,1,-1,-1));
    QPen pen(Qt::darkGreen);
    pen.setWidth(2);
    painter->setPen(pen);

    painter->setRenderHint(QPainter::Antialiasing);

    for (int i = 0; i < v_isolines.size(); ++i){
        for (int j = 0; j < v_isolines[i].size(); ++j){
            double dx1 = v_isolines[i][j].first.x() - settings.minX;
            double dy1 = v_isolines[i][j].first.y() - settings.minY;

            double x1 = rect.left() + (dx1*(rect.width()-1)/settings.spanX());
            double y1 = rect.bottom() - (dy1*(rect.height()-1)/settings.spanY());


            double dx2 = v_isolines[i][j].second.x() - settings.minX;
            double dy2 = v_isolines[i][j].second.y() - settings.minY;

            double x2 = rect.left() + (dx2*(rect.width()-1)/settings.spanX());
            double y2 = rect.bottom() - (dy2*(rect.height()-1)/settings.spanY());

            painter->drawLine(QPointF(x1,y1),QPointF(x2,y2));
        }
    }
}

void contourWidget::addIsoline(const QVector<QPair<QPointF,QPointF>> &isoline)
{
    if (!isoline.size()) return;

    v_isolines.append(isoline);
    update();
}

void contourWidget::addIsolines(const QVector<QVector<QPair<QPointF, QPointF> > > &isolines)
{
    if (!isolines.size()) return;

    bool updateNedded = false;

    for (int i = 0; i < isolines.size(); ++i){
        if (isolines[i].size()){
           v_isolines.append(isolines[i]);
           updateNedded = true;
        }
    }

    if (updateNedded) update();
}

void contourWidget::clearIsolines()
{
    if (v_isolines.size()){
        v_isolines.clear();
        update();
    }
}
PlotSettings::PlotSettings()
{
    minX = 0.0;
    maxX = 10.0;
    numXTicks = 5;

    minY = 0.0;
    maxY = 10.0;
    numYTicks = 5;
}

void PlotSettings::scroll(int dx, int dy)
{
    double stepX = spanX() / numXTicks;
    minX += dx * stepX;
    maxX += dx * stepX;

    double stepY = spanY() / numYTicks;
    minY += dy * stepY;
    maxY += dy * stepY;

}

void PlotSettings::adjust()
{
    adjustAxis(minX, maxX, numXTicks);
    adjustAxis(minY, maxY, numYTicks);
}

void PlotSettings::adjustAxis(double &min, double &max, int &numTicks)
{

    const int MinTicks = 4;
    double grossStep = (max - min) / MinTicks;
    double step = std::pow(10.0, std::floor(std::log10(grossStep)));

    if (5 * step < grossStep) {
        step *= 5;
    } else if (2 * step < grossStep) {
        step *= 2;
    }
    numTicks = int(std::ceil(max / step) - std::floor(min / step));
    if (numTicks < MinTicks)
        numTicks = MinTicks;
    //min = std::ceil(min / step) * step;
    //max = std::floor(max / step) * step;

}
