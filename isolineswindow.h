#ifndef ISOLINESWINDOW_H
#define ISOLINESWINDOW_H

#include <QMainWindow>
#include "contourwidget.h"

class QComboBox;
class QPushButton;
class QTextEdit;
class QGroupBox;
class QDoubleSpinBox;
class QLabel;
class isolinesLineEdit;

class isolinesWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit isolinesWindow(QWidget *parent = nullptr);

    template<typename T>
    void computeIsolines(T &func);

signals:

public slots:
   void graficarNuevasIsolineas();
   void agregarNuevasIsolineas();
private:
    contourWidget *m_contourWidget;
    QTextEdit *m_resultsDisplay;

    QComboBox *m_optTypeComboBox;

    QGroupBox *m_withOutConstGroupBox;
    QComboBox *m_problemWORComboBox;
    QComboBox *m_methodWORComboBox;

    QGroupBox *m_withConstGroupBox;
    QComboBox *m_problemWRComboBox;
    QComboBox *m_methodWRComboBox;
    //QAction *m_qp_eqConstAction;
    //QAction *m_qp_actSetAction;
    //QAction *m_qp_gradProjAction;
    //QAction *m_interiorPointAction;

    QGroupBox *m_isolinesGroupBox;
    QLabel *m_XLabel;
    QLabel *m_YLabel;
    QLabel *m_rangoLabel;
    QLabel *m_stepLabel;
    QLabel *m_isolinesLabel;

    QDoubleSpinBox *m_xbSpinButton;
    QDoubleSpinBox *m_xeSpinButton;
    QDoubleSpinBox *m_ybSpinButton;
    QDoubleSpinBox *m_yeSpinButton;
    QDoubleSpinBox *m_stepSpinBox;

    isolinesLineEdit *m_isolinesLineEdit;

    QAction *m_obtenerIsolineasAction;
    QAction *m_agregarIsolineasAction;
    QMenu *m_isolineasMenu;
    QPushButton *m_getIsolinesButton;

    QPushButton *m_minimizePushButton;

    QVector<double> m_isovalues;
};

#endif // ISOLINESWINDOW_H
