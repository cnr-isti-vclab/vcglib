#ifndef TRACKBALL_H
#define TRACKBALL_H

#include <vcg/math/similar.h>

namespace vcg {

class Trackball {
public:
  Similarf track;
  Similarf local;

  Trackball();
  void SetIdentity();
  void SetPosition(const Similarf &local, int millisec = 0);

  //operating
  void GetView();
  void Apply();
  void Draw();
  void Reset();

  //interface
  void MouseDown(int x, int y, Button button);
	void MouseMove(int x, int y); 
	void MouseUp(int x, int y, Button button); 
	void MouseWheel(Button notch);
  void ButtonUp(Button button);
  void ButtonDown(Button button);

  //spinning interface
  void SetSpinnable(bool on);
  bool IsSpinnable();  
  bool SetSpinning(Quaternionf &spin);
  void StopSpinning();
  bool IsSpinning();  

  //interfaccia navigation:
  void Back();
  void Forward();
  void Home();
  void HistorySize(int lenght);

  enum { LOCAL, VIEW, SCREEN };
  
  enum {
		ROTATE = 0, ROTATE_G = 1,                             //really makes sense only in VIEW system
		ROTATE_X = 2, ROTATE_Y = 3, ROTATE_Z = 4,		    // Axis Constrained Rotation 
		DRAG_X = 5,   DRAG_Y = 6,   DRAG_Z = 7,					// Drag constrained to an axis (trackball axis)
		DRAG_XY = 8,  DRAG_YZ = 9,  DRAG_XZ = 10,				// Drag constrained to a plane
		SCALE = 11,                                  //scale respect to center of trackball
		NONE = 12                                    //disable trackball
  };
	enum Button { BUTTON_LEFT = 1, BUTTON_MIDDLE = 2, BUTTON_RIGHT = 4, WHEEL = 8,
		            KEY_SHIFT = 16, KEY_CTRL = 32, KEY_ALT = 64, HANDLE = 128 };
protected:
  Camera camera;

  TrackMode *current;

  Similarf &last();
  int last_x, last_y;
  bool dragging;
  int button_mask;

  Quaternionf spin;
  bool spinnable;
  bool spinning;

  std::list<Similarf> history;
};

}//namespace

#endif