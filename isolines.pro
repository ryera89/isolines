#-------------------------------------------------
#
# Project created by QtCreator 2017-12-21T16:22:44
#
#-------------------------------------------------

QT       += core gui charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = isolines
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS
CONFIG += c++11

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += \
        main.cpp \
    linalg/matrix_ref1d.cpp \
    linalg/matrix_ref2d.cpp \
    numopt/pointsender.cpp \
    contourwidget.cpp \
    isolineswindow.cpp \
    isolineslineedit.cpp

HEADERS += \
    examples/ejem_16_04_interior_point.h \
    examples/ejemplo_juanma.h \
    examples/ejemplo_rosenbrock_interior_point.h \
    examples/example_convex_qp_active_set.h \
    examples/functors.h \
    examples/interior_point_example.h \
    examples/linesearchnumtestfunc.h \
    examples/prob_16_11_interior_point.h \
    linalg/cholesky.h \
    linalg/directsolvers.h \
    linalg/gauss2x2.h \
    linalg/givens_rotations.h \
    linalg/householdervec.h \
    linalg/itermethods.h \
    linalg/lbldcmp.h \
    linalg/ludcmp.h \
    linalg/mathmatrix.h \
    linalg/mathvector.h \
    linalg/matrix_ref1d.h \
    linalg/matrix_ref2d.h \
    linalg/my_concepts.h \
    linalg/qr.h \
    linalg/utils.h \
    numopt/cgfrpr.h \
    numopt/cghz.h \
    numopt/convanalizer.h \
    numopt/interiorpoint.h \
    numopt/interpolmin.h \
    numopt/linesearch.h \
    numopt/numopt.h \
    numopt/qnbfgs.h \
    numopt/qp_gpm_bc.h \
    numopt/qp_lin_convex_act_set.h \
    numopt/qp_lin_eqc.h \
    numopt/steepestdescent.h \
    numopt/trustregion.h \
    numopt/pointsender.h \
    examples/isolines_functions.h \
    marching_squares.h \
    interp_1d.h \
    contourwidget.h \
    isolineswindow.h \
    isolineslineedit.h

FORMS +=
