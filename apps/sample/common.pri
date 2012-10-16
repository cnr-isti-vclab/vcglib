DEPENDPATH += . ../../..
INCLUDEPATH += . ../../..
CONFIG += console stl
TEMPLATE = app
# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
