#include "isolineslineedit.h"
#include <QKeyEvent>

isolinesLineEdit::isolinesLineEdit(QWidget *parent):QLineEdit(parent)
{

}

QVector<double> isolinesLineEdit::getIsolinesValues()
{
    QString stringValues = this->text();

    QStringList values = stringValues.split(";");
    QVector<double> isovalues;
    isovalues.reserve(values.size());
    for (int i = 0; i < values.size(); ++i){
        bool ok;
        double val = values.at(i).toDouble(&ok);
        if (ok) isovalues.append(val);
    }

    return isovalues;
}

void isolinesLineEdit::keyPressEvent(QKeyEvent *event)
{
    int key = event->key();

    if (key != Qt::Key_0 && key != Qt::Key_1
            && key != Qt::Key_2 && key != Qt::Key_3
            && key != Qt::Key_4 && key != Qt::Key_5
            && key != Qt::Key_6 && key != Qt::Key_7
            && key != Qt::Key_8 && key != Qt::Key_9
            && key != Qt::Key_Semicolon && key != Qt::Key_Enter
            && key != Qt::Key_Period
            && key != Qt::Key_Minus && key != Qt::Key_Backspace
            && key != Qt::Key_Delete){
        return;
    }
    if ((prevKey == Qt::Key_Period && key == Qt::Key_Semicolon) ||
            (prevKey == Qt::Key_Semicolon && key == Qt::Key_Semicolon) ||
            (prevKey == Qt::Key_Period && key == Qt::Key_Period) ||
            (prevKey == Qt::Key_Semicolon && key == Qt::Key_Period) ||
            (prevKey == Qt::Key_Period && key == Qt::Key_Minus) ||
            (prevKey == Qt::Key_Minus && key == Qt::Key_Semicolon)){
        return;
    }
    if ((prevKey == Qt::Key_0 || prevKey == Qt::Key_1
            || prevKey == Qt::Key_2 || prevKey == Qt::Key_3
            || prevKey == Qt::Key_4 || prevKey == Qt::Key_5
            || prevKey == Qt::Key_6 || prevKey == Qt::Key_7
            || prevKey == Qt::Key_8 || prevKey == Qt::Key_9) && key == Qt::Key_Minus){
        return;
    }
    prevKey = key;
    QLineEdit::keyPressEvent(event);
}
