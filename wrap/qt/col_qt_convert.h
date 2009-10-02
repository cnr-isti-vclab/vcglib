#ifndef COL_QT_CONVERT_H_
#define COL_QT_CONVERT_H_

#include <QColor>
#include "../../vcg/space/color4.h"

namespace vcg
{
	class ColorConverter
	{
	public:
		template<typename T>
		static vcg::Color4<T> convertQColorToColor4(const QColor& col)
		{
			return vcg::Color4<T>(T(col.red()),T(col.green()),T(col.blue()),T(col.alpha()));
		}

		template<typename T>
		static QColor convertColor4ToQColor(const vcg::Color4<T>& col)
		{
			return QColor(int(col.X()),int(col.Y()),int(col.Z()),int(col.W()));
		}
	};
}

#endif