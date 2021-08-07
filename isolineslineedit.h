#ifndef ISOLINESLINEEDIT_H
#define ISOLINESLINEEDIT_H

#include <QLineEdit>

class isolinesLineEdit : public QLineEdit
{
    Q_OBJECT
public:
    isolinesLineEdit(QWidget *parent = nullptr);

    QVector<double> getIsolinesValues();

protected:
    void keyPressEvent(QKeyEvent *event);

private:
    int prevKey;
};

#endif // ISOLINESLINEEDIT_H
