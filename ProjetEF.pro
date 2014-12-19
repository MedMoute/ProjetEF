#-------------------------------------------------
#
# Project created by QtCreator 2014-10-25T20:07:58
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = ProjetEF

CONFIG   += console
CONFIG   += app_bundle


TEMPLATE = app


SOURCES += \
    src/maillage.cpp \
    src/probleme.cpp \
    src/main.cpp \
    src/parallel.cpp \
    src/nonParallel.cpp

HEADERS += \
    include/probleme.h \
    include/maillage.h \
    include/parallel.h \
    include/nonParallel.h \
    include/param.h

OTHER_FILES += \
    fichierTest/testpart.msh

#Add path to Eigen
INCLUDEPATH += ./eigen
INCLUDEPATH += "/Users/pierrefournier/opt/openmpi/include"
LIBS += "/Users/pierrefournier/opt/openmpi/lib"
