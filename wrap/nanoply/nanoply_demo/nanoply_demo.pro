DEPENDPATH += ..
INCLUDEPATH += ..
CONFIG += console stl c++11
CONFIG -= qt
TEMPLATE = app
# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
SOURCES += main.cpp
