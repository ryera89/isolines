#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "isolines_widget.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    isolines_widget* rosenbrockIsoView(){return rosenbrockIso;}
    isolines_widget* ejemplo_1604IsoView(){return ejemplo_1604Iso;}
    isolines_widget* ejemplo_1611IsoView(){return ejemplo_1611Iso;}
    isolines_widget* ejemplo_juanmaIsoView(){return ejemplo_juanmaIso;}

private:
    Ui::MainWindow *ui;
    isolines_widget *rosenbrockIso;
    isolines_widget *ejemplo_1604Iso;
    isolines_widget *ejemplo_1611Iso;
    isolines_widget *ejemplo_juanmaIso;

};

#endif // MAINWINDOW_H
