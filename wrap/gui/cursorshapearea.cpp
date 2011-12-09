#include "cursorshapearea.h"

CursorShapeArea::CursorShapeArea(QDeclarativeItem *parent) :
  QDeclarativeItem(parent)
{
}

Qt::CursorShape CursorShapeArea::cursorShape() const
{
  return cursor().shape();
}

void CursorShapeArea::setCursorShape(Qt::CursorShape cursorShape)
{
  if (cursor().shape() == cursorShape)
    return;

  setCursor(cursorShape);
  emit cursorShapeChanged();
}
