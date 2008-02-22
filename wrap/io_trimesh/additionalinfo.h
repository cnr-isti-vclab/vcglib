#ifndef __VCGLIB_ADDITIONAL_INFO
#define __VCGLIB_ADDITIONAL_INFO

class AdditionalInfo
{
protected:
	AdditionalInfo()
	{
	}
public:
	unsigned int numvert;
	unsigned int numface;
	int mask;

	virtual ~AdditionalInfo()
	{
	}
};

#endif
