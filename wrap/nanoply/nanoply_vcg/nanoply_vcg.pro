INCLUDEPATH += ../../../ ../../../eigenlib
CONFIG += console stl c++17
CONFIG += qt
TEMPLATE = app
SOURCES += main.cpp
QT	    += core

# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
