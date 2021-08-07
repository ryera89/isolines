#ifndef POINTSENDER_H
#define POINTSENDER_H

#include <QObject>

class pointSender : public QObject
{
    Q_OBJECT
public:
    explicit pointSender(QObject *parent = nullptr);

signals:
    void pointUpdated(const QPointF &point);

public slots:
};

#endif // POINTSENDER_H
