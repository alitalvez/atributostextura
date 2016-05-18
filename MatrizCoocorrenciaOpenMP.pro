#-------------------------------------------------
#
# Project created by QtCreator 2015-11-16T14:24:38
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = MatrizCoocorrenciaOpenMP
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    #atCpu.cpp \
    gerarmatriz.cpp \
    hex2int.cxx \
    #matrizcoocorrencia.cpp \
    readimage.cpp \
    haralick.cpp \
    atCpu.cpp \
    matrizcoocorrencia.cpp

HEADERS += \
    #atCpu.h \
    atKernel.cuh \
    gerarmatriz.h \
    globalvar.h \
    hex2int.h \
    #matrizcoocorrencia.h \
    readimage.h \
    tempo.h \
    utilidades.h \
    haralick.h \
    atCpu.h \
    matrizcoocorrencia.h

CONFIG -= qt
#QMAKE_CXXFLAGS_RELEASE = -fopenmp -O3 #-fopt-info-vec-missed
#QMAKE_LFLAGS_RELEASE = -fopenmp -O3
QMAKE_CXXFLAGS+= -fopenmp
#QMAKE_LFLAGS +=  -openmp
LIBS += -fopenmp

DESTDIR     = $$system(pwd)
OBJECTS_DIR = $$DESTDIR/Obj

DISTFILES += \
    D4500.LEFT_CC.raw \
    gmon.out \
    analise.txt \
    Makefile
