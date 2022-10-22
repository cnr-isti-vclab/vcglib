#ifndef __VCGLIB_COLORMAP_H
#define __VCGLIB_COLORMAP_H

#include <vcg/space/color4.h>
#include <map>
#include <vector>

// Colormaps sampled from Matplotlib
// https://matplotlib.org

namespace vcg {

enum ColorMap {
	RGB = 0,
	Viridis = 1,
	Plasma = 2,
	Cividis = 3,
	Turbo = 4,
	RdPu = 5,
};

static std::vector<ColorMap> ColorMapEnums = {
	RGB,
	Viridis,
	Plasma,
	Cividis,
	Turbo,
	RdPu,
};

static std::map<ColorMap, std::vector<vcg::Color4b>> colorMaps = {
	{
		ColorMap::Viridis,
		{
			vcg::Color4b(68, 1, 84, 255),
			vcg::Color4b(70, 12, 95, 255),
			vcg::Color4b(71, 24, 106, 255),
			vcg::Color4b(72, 34, 115, 255),
			vcg::Color4b(70, 45, 124, 255),
			vcg::Color4b(68, 55, 129, 255),
			vcg::Color4b(65, 65, 134, 255),
			vcg::Color4b(61, 74, 137, 255),
			vcg::Color4b(57, 84, 139, 255),
			vcg::Color4b(53, 92, 140, 255),
			vcg::Color4b(49, 100, 141, 255),
			vcg::Color4b(46, 108, 142, 255),
			vcg::Color4b(42, 117, 142, 255),
			vcg::Color4b(39, 124, 142, 255),
			vcg::Color4b(36, 132, 141, 255),
			vcg::Color4b(34, 139, 141, 255),
			vcg::Color4b(31, 148, 139, 255),
			vcg::Color4b(30, 155, 137, 255),
			vcg::Color4b(31, 163, 134, 255),
			vcg::Color4b(36, 170, 130, 255),
			vcg::Color4b(46, 178, 124, 255),
			vcg::Color4b(57, 185, 118, 255),
			vcg::Color4b(71, 192, 110, 255),
			vcg::Color4b(87, 198, 101, 255),
			vcg::Color4b(107, 205, 89, 255),
			vcg::Color4b(126, 210, 78, 255),
			vcg::Color4b(146, 215, 65, 255),
			vcg::Color4b(167, 219, 51, 255),
			vcg::Color4b(191, 223, 36, 255),
			vcg::Color4b(212, 225, 26, 255),
			vcg::Color4b(233, 228, 25, 255),
			vcg::Color4b(253, 231, 36, 255),
		}
	},
	{
		ColorMap::Plasma,
		{
			vcg::Color4b(12, 7, 134, 255),
			vcg::Color4b(33, 5, 143, 255),
			vcg::Color4b(49, 4, 150, 255),
			vcg::Color4b(63, 3, 156, 255),
			vcg::Color4b(78, 2, 161, 255),
			vcg::Color4b(90, 0, 165, 255),
			vcg::Color4b(103, 0, 167, 255),
			vcg::Color4b(115, 0, 168, 255),
			vcg::Color4b(129, 4, 167, 255),
			vcg::Color4b(140, 10, 164, 255),
			vcg::Color4b(151, 19, 160, 255),
			vcg::Color4b(162, 28, 154, 255),
			vcg::Color4b(173, 38, 146, 255),
			vcg::Color4b(182, 47, 139, 255),
			vcg::Color4b(190, 56, 131, 255),
			vcg::Color4b(198, 65, 124, 255),
			vcg::Color4b(207, 75, 116, 255),
			vcg::Color4b(214, 85, 109, 255),
			vcg::Color4b(220, 94, 102, 255),
			vcg::Color4b(227, 103, 95, 255),
			vcg::Color4b(233, 114, 87, 255),
			vcg::Color4b(238, 124, 80, 255),
			vcg::Color4b(243, 134, 73, 255),
			vcg::Color4b(246, 145, 66, 255),
			vcg::Color4b(250, 157, 58, 255),
			vcg::Color4b(252, 169, 52, 255),
			vcg::Color4b(253, 181, 45, 255),
			vcg::Color4b(253, 193, 40, 255),
			vcg::Color4b(251, 208, 36, 255),
			vcg::Color4b(248, 221, 36, 255),
			vcg::Color4b(244, 234, 38, 255),
			vcg::Color4b(239, 248, 33, 255),
		}
	},
	{
		ColorMap::Cividis,
		{
			vcg::Color4b(0, 34, 77, 255),
			vcg::Color4b(0, 40, 91, 255),
			vcg::Color4b(0, 45, 105, 255),
			vcg::Color4b(4, 50, 112, 255),
			vcg::Color4b(28, 56, 110, 255),
			vcg::Color4b(40, 62, 109, 255),
			vcg::Color4b(50, 68, 108, 255),
			vcg::Color4b(59, 73, 107, 255),
			vcg::Color4b(69, 79, 107, 255),
			vcg::Color4b(77, 85, 108, 255),
			vcg::Color4b(84, 90, 108, 255),
			vcg::Color4b(91, 96, 110, 255),
			vcg::Color4b(99, 102, 111, 255),
			vcg::Color4b(106, 108, 113, 255),
			vcg::Color4b(113, 114, 115, 255),
			vcg::Color4b(120, 120, 118, 255),
			vcg::Color4b(128, 126, 120, 255),
			vcg::Color4b(135, 132, 120, 255),
			vcg::Color4b(143, 138, 119, 255),
			vcg::Color4b(151, 144, 118, 255),
			vcg::Color4b(160, 151, 117, 255),
			vcg::Color4b(168, 158, 115, 255),
			vcg::Color4b(176, 164, 112, 255),
			vcg::Color4b(184, 171, 109, 255),
			vcg::Color4b(194, 178, 105, 255),
			vcg::Color4b(202, 185, 100, 255),
			vcg::Color4b(211, 192, 95, 255),
			vcg::Color4b(219, 199, 89, 255),
			vcg::Color4b(229, 207, 80, 255),
			vcg::Color4b(238, 215, 71, 255),
			vcg::Color4b(248, 222, 59, 255),
			vcg::Color4b(253, 231, 55, 255),
		}
	},
	{
		ColorMap::Turbo,
		{
			vcg::Color4b(48, 18, 59, 255),
			vcg::Color4b(57, 41, 114, 255),
			vcg::Color4b(64, 64, 161, 255),
			vcg::Color4b(68, 86, 199, 255),
			vcg::Color4b(70, 109, 230, 255),
			vcg::Color4b(70, 130, 248, 255),
			vcg::Color4b(64, 150, 254, 255),
			vcg::Color4b(52, 170, 248, 255),
			vcg::Color4b(37, 192, 230, 255),
			vcg::Color4b(26, 209, 210, 255),
			vcg::Color4b(24, 224, 189, 255),
			vcg::Color4b(34, 235, 169, 255),
			vcg::Color4b(59, 244, 141, 255),
			vcg::Color4b(89, 251, 114, 255),
			vcg::Color4b(120, 254, 89, 255),
			vcg::Color4b(149, 254, 68, 255),
			vcg::Color4b(174, 249, 55, 255),
			vcg::Color4b(195, 241, 51, 255),
			vcg::Color4b(214, 229, 53, 255),
			vcg::Color4b(231, 215, 56, 255),
			vcg::Color4b(244, 196, 58, 255),
			vcg::Color4b(251, 179, 54, 255),
			vcg::Color4b(254, 158, 46, 255),
			vcg::Color4b(252, 134, 36, 255),
			vcg::Color4b(246, 107, 24, 255),
			vcg::Color4b(237, 85, 15, 255),
			vcg::Color4b(226, 66, 9, 255),
			vcg::Color4b(212, 50, 5, 255),
			vcg::Color4b(192, 35, 2, 255),
			vcg::Color4b(172, 22, 1, 255),
			vcg::Color4b(148, 12, 1, 255),
			vcg::Color4b(122, 4, 2, 255),
		}
	},
	{
		ColorMap::RdPu,
		{
			vcg::Color4b(255, 247, 243, 255),
			vcg::Color4b(254, 241, 237, 255),
			vcg::Color4b(253, 235, 231, 255),
			vcg::Color4b(253, 229, 226, 255),
			vcg::Color4b(252, 223, 219, 255),
			vcg::Color4b(252, 216, 212, 255),
			vcg::Color4b(252, 209, 205, 255),
			vcg::Color4b(252, 202, 198, 255),
			vcg::Color4b(251, 194, 191, 255),
			vcg::Color4b(251, 184, 188, 255),
			vcg::Color4b(250, 175, 185, 255),
			vcg::Color4b(250, 165, 182, 255),
			vcg::Color4b(249, 153, 178, 255),
			vcg::Color4b(248, 139, 173, 255),
			vcg::Color4b(248, 125, 168, 255),
			vcg::Color4b(247, 111, 163, 255),
			vcg::Color4b(243, 96, 159, 255),
			vcg::Color4b(236, 83, 157, 255),
			vcg::Color4b(230, 70, 154, 255),
			vcg::Color4b(223, 57, 152, 255),
			vcg::Color4b(212, 42, 146, 255),
			vcg::Color4b(200, 30, 140, 255),
			vcg::Color4b(189, 17, 134, 255),
			vcg::Color4b(177, 4, 127, 255),
			vcg::Color4b(162, 1, 124, 255),
			vcg::Color4b(149, 1, 122, 255),
			vcg::Color4b(136, 1, 121, 255),
			vcg::Color4b(123, 1, 119, 255),
			vcg::Color4b(109, 0, 115, 255),
			vcg::Color4b(97, 0, 112, 255),
			vcg::Color4b(85, 0, 109, 255),
			vcg::Color4b(73, 0, 106, 255),
		}
	},
};

inline vcg::Color4b GetColorMapping(double v, double minv, double maxv, ColorMap cmap)
{
	if (cmap == ColorMap::RGB) {
		vcg::Color4b c;
		c.SetColorRamp(float(minv), float(maxv), float(v));
		return c;
	}

	int sz = int(colorMaps[cmap].size());

	if (v > maxv)
		v = maxv;
	if (v < minv)
		v = minv;

	v = (v - minv) / (maxv - minv);

	double v0 = v * sz;
	int n0 = int(v0);

	if (n0 < 0)
		return colorMaps[cmap].front();
	if (n0 >= sz - 1)
		return colorMaps[cmap].back();

	int n1 = n0 + 1;

	double fract = v0 - n0;

	vcg::Color4b c0 = colorMaps[cmap][n0];
	vcg::Color4b c1 = colorMaps[cmap][n1];

	c0.lerp(c0, c1, fract);

	return c0;
}

} // namespace vcg

#endif
