DEPENDPATH += \
	. \
	../../..

INCLUDEPATH += \
	. \
	../../.. \
	../../../eigenlib


CONFIG += c++11
TEMPLATE = app

# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle

QMAKE_CXXFLAGS += -std=c++11

unix {
	CONFIG(release, debug|release) {
		DEFINES *= NDEBUG
	}
}
