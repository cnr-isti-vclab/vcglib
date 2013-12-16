#ifndef DEVICE_TO_LOGICAL_H
#define DEVICE_TO_LOGICAL_H
#include <QWidget>
#include <QPainter>
template < class ValueType>
inline ValueType QTLogicalToDevice( QWidget *qw, const ValueType &value)
{
#if QT_VERSION >= 0x050000
  return value*qw->devicePixelRatio() ;
#else
  Q_UNUSED(qw);
  return value;
#endif
}

template < class ValueType>
inline ValueType QTLogicalToDevice( QPainter *qp, const ValueType &value)
{
#if QT_VERSION >= 0x050000
  return value*qp->device()->devicePixelRatio() ;
#else
  Q_UNUSED(qp);
  return value;
#endif
}

template < class ValueType>
inline ValueType QTDeviceToLogical( QPainter *qp, const ValueType &value)
{
#if QT_VERSION >= 0x050000
  return value/qp->device()->devicePixelRatio() ;
#else
  Q_UNUSED(qp);
  return value;
#endif
}

#endif // DEVICE_TO_LOGICAL_H
