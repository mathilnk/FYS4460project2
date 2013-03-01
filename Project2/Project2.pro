TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    atom.cpp \
    cell.cpp \
    cellcontainer.cpp \
    cellsolver.cpp \
    thermostattest.cpp\
    normal.cpp

HEADERS += \
    atom.h \
    cell.h \
    cellcontainer.h \
    cellsolver.h \
    thermostattest.h \
    normal.hpp

LIBS +=-fopenmp

COMMON_CXXFLAGS +=  -fopenmp
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS
