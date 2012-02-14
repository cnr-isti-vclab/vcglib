#ifndef ANTTWEAKBARMAPPER_H
#define ANTTWEAKBARMAPPER_H

void TW_CALL CopyCDStringToClient(char **destPtr, const char *src);
void TW_CALL CopyStdStringToClient(std::string& destClientString, const std::string& srcLibraryString);

TwMouseButtonID Qt2TwMouseButtonId(QMouseEvent *e);
int TwMousePressQt(QMouseEvent *e);
int TwMouseReleaseQt(QMouseEvent *e);
int TwKeyPressQt(QKeyEvent *e);

#endif // ANTTWEAKBARMAPPER_H
