/****************************************************************************
* NanoPLY                                                                   *
* NanoPLY is a C++11 header-only library to read and write PLY file         *
*                                                                           *
* Copyright(C) 2014                                                         *
* Visual Computing Lab                                                      *
* ISTI - Italian National Research Council                                  *
*                                                                           *
* This Source Code Form is subject to the terms of the Mozilla Public       *
* License, v. 2.0. If a copy of the MPL was not distributed with this       *
* file, You can obtain one at http://mozilla.org/MPL/2.0/.                  *
*                                                                           *
****************************************************************************/


#ifndef NANOPLY_HPP
#define NANOPLY_HPP

#include <vector>
#include <tuple>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <cstdio>
#include <cmath>
#include <limits>
#include <stdint.h>

// Avoid conflicting declaration of min/max macros in windows headers
#if !defined(NOMINMAX) && (defined(_WIN32) || defined(_WIN32_)  || defined(WIN32) || defined(_WIN64))
# define NOMINMAX
# ifdef max
#  undef   max
#  undef   min
# endif
#endif

#define USE_H4D_OUTPUT

namespace nanoply
{

	/**
	 * @cond HIDDEN_SYMBOLS
	 */
    template < size_t T> struct SizeT {};
	/**
	 * @endcond
	 */

	/** Error Type.
	 *	Error type returned by the open of a PLY file.
	 */
    typedef enum NNP_ERROR
    {
        NNP_OK                = 0x0000,		/**< No error. */
        NNP_UNABLE_TO_OPEN    = 0x0001,		/**< The file cannot be opend. */
        NNP_MISSING_HEADER    = 0x0002,		/**< The file does not contain a valid PLY header. */
        NNP_MISSING_FORMAT    = 0x0004		/**< The file has an invalid internal format. */
    } ErrorCode;


	/** PLY Entity.
	 *  Property that can be saved in a PLY file.
	 */
    typedef enum NNP_ENTITY
    {
        NNP_UNKNOWN_ENTITY = 0x0000,	/**< Unknown property. */
        NNP_PX             = 0x0001,	/**< Position x cordinate. */
        NNP_PY             = 0x0002,	/**< Position y cordinate. */
        NNP_PZ             = 0x0004, 	/**< Position z cordinate. */
        NNP_PXYZ           = 0x0007, 	/**< Position (x, y, z). */

        NNP_NX             = 0x0010, 	/**< Normal x component. */
        NNP_NY             = 0x0020, 	/**< Normal y component. */
        NNP_NZ             = 0x0040, 	/**< Normal z component. */
        NNP_NXYZ           = 0x0070,  	/**< Normal (x, y, z). */

        NNP_CR             = 0x0100,  	/**< Color red component. */
        NNP_CG             = 0x0200,  	/**< Color green component. */
        NNP_CB             = 0x0400,  	/**< Color blue component. */
        NNP_CRGB           = 0x0700,  	/**< Color RGB. */
        NNP_CA             = 0x0800,  	/**< Color alpha component. */
        NNP_CRGBA          = 0x0F00,  	/**< Color RGBA. */

        NNP_DENSITY        = 0x1000,   	/**< Density or Radius property. */
        NNP_SCALE          = 0x2000,  	/**< Scale property. */
        NNP_QUALITY		   = 0x4000,  	/**< Quality property. */
        NNP_REFLECTANCE	   = 0x8000,  	/**< Refelctance property. */
		NNP_BITFLAG		   = 0x10000,   /**< Bit flags. */

        NNP_VERTEX_LIST    = 0x20000  	/**< List of vertec index. */
    } PlyEntity;


	/** PLY Type.
	 *  Type of a PLY property.
	 */
    typedef enum NNP_PLYTYPE
    {
        NNP_UNKNOWN_TYPE      = 0x0000, /**< Unknown type. */
        NNP_FLOAT32           = 0x0001, /**< Float. */
        NNP_FLOAT64           = 0x0002, /**< Double. */
        NNP_INT8              = 0x0010, /**< Char. */
        NNP_INT16             = 0x0020, /**< Short. */
        NNP_INT32             = 0x0040, /**< Int. */
        NNP_UINT8             = 0x0100, /**< Unsigned Char. */
        NNP_UINT16            = 0x0200, /**< Unsigned Short. */
        NNP_UINT32            = 0x0400, /**< Unsigned Int. */
        NNP_LIST_UINT8_UINT32 = 0x1000, /**< List (size Unsigned Char) of Unsigned Int.  */
        NNP_LIST_INT8_UINT32  = 0x2000, /**< List (size Char) of Unsigned Int.  */
        NNP_LIST_UINT8_INT32  = 0x4000, /**< List (size Unsigned Char) of Int.  */
        NNP_LIST_INT8_INT32   = 0x8000  /**< List (size Char) of Int. */
    } PlyType;


	/** PLY Property.
	 *  Define a PLY property (entity and type).
	 */
    class PlyProperty
    {
    public:
        PlyType type;		/**< Property type. */
        PlyEntity elem;		/**< Property entity. */


		/** 
		 * Constructor that sets the type and the entity of the PLY property.
		 *
		 * @param _t Property type.
		 * @param _e Property entity.
		 */
        PlyProperty(PlyType _t, PlyEntity _e):type(_t),elem(_e){}

		/** 
		 * Get the description string of the property entity. 
		 *
		 * @return Description string of the property entity. 
		 */
        const char* EntityStr();

		/**
		* Get the name of the property entity.
		*
		* @return Name of the property entity.
		*/
		const char* EntityName();

		/** 
		 * Get the description string of the property type. 
		 *
		 * @return Description string of the property type. 
		 */
		const char* TypeStr();

		/** 
		 * Get the size in byte of the property type. 
		 *
		 * @return Size in byte of the property type. 
		 */
        int TypeSize();

		/** 
		 * Skip the property in an Ascii file. 
		 *
		 * @param *fp	Pointer to the opened file. 
		 * @return		If successful returns true. Otherwise, it returns false. 
		 */
        bool SkipAsciiPropertyInFile(FILE *fp);
        
		/** 
		 * Skip the property in a binary file. 
		 *
		 * @param *fp	Pointer to the opened file. 
		 * @return		If successful returns true. Otherwise, it returns false. 
		 */
		bool SkipBinaryPropertyInFile(FILE *fp);

		/** 
		 * Write the property in the header of the PLY file. 
		 *
		 * @param *fp	Pointer to the file. 
		 * @return		If successful returns true. Otherwise, it returns false. 
		 */
        bool WriteHeader(FILE *fp);
    };

    const char* PlyProperty::EntityStr()
    {
        switch (this->elem)
        {
        case NNP_UNKNOWN_ENTITY  :    return "NNP_UNKNOWN_ENTITY";
        case NNP_PX              :    return "NNP_PX            ";
        case NNP_PY              :    return "NNP_PY            ";
        case NNP_PZ              :    return "NNP_PZ            ";
        case NNP_PXYZ            :    return "NNP_PXYZ          ";

        case NNP_NX              :    return "NNP_NX            ";
        case NNP_NY              :    return "NNP_NY            ";
        case NNP_NZ              :    return "NNP_NZ            ";
        case NNP_NXYZ            :    return "NNP_NXYZ          ";

        case NNP_CR              :    return "NNP_CR            ";
        case NNP_CG              :    return "NNP_CG            ";
        case NNP_CB              :    return "NNP_CB            ";
        case NNP_CRGB            :    return "NNP_CRGB          ";
        case NNP_CA              :    return "NNP_CA            ";
        case NNP_CRGBA			 :    return "NNP_CRGBA         ";

        case NNP_DENSITY         :    return "NNP_DENSITY       ";
        case NNP_SCALE           :    return "NNP_SCALE         ";
        case NNP_QUALITY		 :    return "NNP_QUALITY       ";
        case NNP_REFLECTANCE	 :    return "NNP_REFLECTANCE   ";

        case NNP_VERTEX_LIST	 :    return "NNP_VERTEX_LIST   ";

        default: assert(0);
            break;
        }
        return 0;
    }


	const char* PlyProperty::EntityName()
	{
#ifdef USE_H4D_OUTPUT
		switch (this->elem)
		{
		case NNP_UNKNOWN_ENTITY	:    return "unknown";
		case NNP_PX				:    return "x";
		case NNP_PY				:    return "y";
		case NNP_PZ				:    return "z";
		case NNP_PXYZ			:    return "x y z";

		case NNP_NX				:    return "nx";
		case NNP_NY				:    return "ny";
		case NNP_NZ				:    return "nz";
		case NNP_NXYZ			:    return "nx ny nz";

		case NNP_CR				:    return "red";
		case NNP_CG				:    return "green";
		case NNP_CB				:    return "blue";
		case NNP_CRGB			:    return "rgb";
		case NNP_CA				:    return "alpha";
		case NNP_CRGBA			:    return "rgba";

		case NNP_DENSITY		:    return "radius";
		case NNP_SCALE			:    return "scale";
		case NNP_QUALITY		:    return "quality";
		case NNP_REFLECTANCE	:    return "reflectance";

		case NNP_VERTEX_LIST	:    return "vertex_indices";

		default: assert(0);
			break;
		}
#else
		switch (this->elem)
		{
		case NNP_UNKNOWN_ENTITY	:    return "unknown";
		case NNP_PX				:    return "x";
		case NNP_PY				:    return "y";
		case NNP_PZ				:    return "z";
		case NNP_PXYZ			:    return "x y z";

		case NNP_NX				:    return "nx";
		case NNP_NY				:    return "ny";
		case NNP_NZ				:    return "nz";
		case NNP_NXYZ			:    return "nx ny nz";

		case NNP_CR				:    return "diffuse_red";
		case NNP_CG				:    return "diffuse_green";
		case NNP_CB				:    return "diffuse_blue";
		case NNP_CRGB			:    return "diffuse_rgb";
		case NNP_CA				:    return "alpha";
		case NNP_CRGBA			:    return "diffuse_rgba";

		case NNP_DENSITY		:    return "radius";
		case NNP_SCALE			:    return "value";
		case NNP_QUALITY		:    return "confidence";
		case NNP_REFLECTANCE	:    return "reflectance";

		case NNP_VERTEX_LIST	:    return "vertex_indices";

		default: assert(0);
			break;
		}
#endif

		return 0;
	}


    const char* PlyProperty::TypeStr()
    {
        switch (this->type)
        {
        case NNP_UNKNOWN_TYPE     :    return "NNP_UNKNOWN_TYPE     ";
        case NNP_FLOAT32          :    return "NNP_FLOAT32          ";
        case NNP_FLOAT64          :    return "NNP_FLOAT64          ";
        case NNP_INT8             :    return "NNP_INT8             ";
        case NNP_INT16            :    return "NNP_INT16            ";
        case NNP_INT32            :    return "NNP_INT32            ";
        case NNP_UINT8            :    return "NNP_UINT8            ";
        case NNP_UINT16           :    return "NNP_UINT16           ";
        case NNP_UINT32           :    return "NNP_UINT32           ";
        case NNP_LIST_UINT8_UINT32:    return "NNP_LIST_UINT8_UINT32";
        case NNP_LIST_INT8_UINT32 :    return "NNP_LIST_INT8_UINT32 ";
        case NNP_LIST_UINT8_INT32 :    return "NNP_LIST_UINT8_INT32 ";
        case NNP_LIST_INT8_INT32  :    return "NNP_LIST_INT8_INT32  ";
        default: assert(0);
            break;
        }
        return 0;
    }

    int PlyProperty::TypeSize()
    {
        switch (this->type)
        {
            case NNP_UNKNOWN_TYPE:
                return 0;
            case NNP_INT8:
            case NNP_UINT8:
                return 1;
            case NNP_INT16:
            case NNP_UINT16:
                return 2;
            case NNP_FLOAT32:
            case NNP_INT32:
            case NNP_UINT32:
                return 4;
            case NNP_FLOAT64:
                return 8;
            case NNP_LIST_UINT8_UINT32:
            case NNP_LIST_INT8_UINT32 :
            case NNP_LIST_UINT8_INT32 :
            case NNP_LIST_INT8_INT32  :
                return 1;
            default: assert(0);
            break;
        }
        return 0;
    }


    bool PlyProperty::SkipAsciiPropertyInFile(FILE *fp)
    {
        int count = 1;
        if (this->elem == NNP_CRGB || this->elem == NNP_NXYZ || this->elem == NNP_PXYZ)
            count = 3;
        else if (this->elem == NNP_CRGBA)
            count = 4;
        switch(type)
        {
            case NNP_INT8:
            case NNP_INT16:
            case NNP_INT32:
            {
                int* temp = new int[count];
                for (int i = 0; i < count ; i++)
                    fscanf(fp, "%d", &temp[i]);
                delete[] temp;
                break;
            }
            case NNP_UINT8:
            case NNP_UINT16:
            case NNP_UINT32:
            {
                unsigned int* temp = new unsigned int[count];
                for (int i = 0; i < count ; i++)
                    fscanf(fp, "%u", &temp[i]);
                delete[] temp;
                break;
            }
            case NNP_FLOAT32:
            {
                float* temp = new float[count];
                for (int i = 0; i < count ; i++)
                    fscanf(fp, "%f", &temp[i]);
                delete[] temp;
                break;
            }
            case NNP_FLOAT64:
            {
                double* temp = new double[count];
                for (int i = 0; i < count ; i++)
                    fscanf(fp, "%f", &temp[i]);
                delete[] temp;
                break;
            }
            case NNP_LIST_UINT8_UINT32:
            case NNP_LIST_INT8_UINT32 :
            {
                fscanf(fp, "%d", &count);
                unsigned int* temp = new unsigned int[count];
                for (int i = 0; i < count ; i++)
                    fscanf(fp, "%u", &temp[i]);
                delete[] temp;
                break;
            }
            case NNP_LIST_UINT8_INT32 :
            case NNP_LIST_INT8_INT32  :
            {
                fscanf(fp, "%d", &count);
                int* temp = new int[count];
                for (int i = 0; i < count ; i++)
                    fscanf(fp, "%d", &temp[i]);
                delete[] temp;
                break;
            }
        }
        return true;
    }

    bool PlyProperty::SkipBinaryPropertyInFile(FILE *fp)
    {
        int count = 1;
        if (this->elem == NNP_CRGB || this->elem == NNP_NXYZ || this->elem == NNP_PXYZ)
            count = 3;
        else if (this->elem == NNP_CRGBA)
            count = 4;
        if (this->type >= NNP_LIST_UINT8_UINT32)
        {
            unsigned char cntList = 0;
            fread(&cntList, sizeof(char), 1, fp);
            fseek(fp, 4 * cntList, SEEK_CUR);
        }
        else
            fseek(fp, this->TypeSize() * count, SEEK_CUR);
        return true;
    }


    bool PlyProperty::WriteHeader(FILE *fp)
    {
        std::string name, type;

        switch (this->type)
        {
            case NNP_UNKNOWN_TYPE     :    type = "unknonw"; break;
            case NNP_FLOAT32          :    type = "float"; break;
            case NNP_FLOAT64          :    type = "double"; break;
            case NNP_INT8             :    type = "char"; break;
            case NNP_INT16            :    type = "short"; break;
            case NNP_INT32            :    type = "int"; break;
            case NNP_UINT8            :    type = "uchar"; break;
            case NNP_UINT16           :    type = "ushort"; break;
            case NNP_UINT32           :    type = "uint"; break;
            case NNP_LIST_UINT8_UINT32:    type = "list uchar uint"; break;
            case NNP_LIST_INT8_UINT32 :    type = "list char uint"; break;
            case NNP_LIST_UINT8_INT32 :    type = "list uchar int"; break;
            case NNP_LIST_INT8_INT32  :    type = "list char int"; break;
        }

		name = this->EntityName();
        switch (this->elem)
        {
            case NNP_PXYZ            :
                {
					name = PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_PX).EntityName();
                    fprintf(fp, "property %s %s\n", type.c_str(), name.c_str());
					name = PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_PY).EntityName();
                    fprintf(fp, "property %s %s\n", type.c_str(), name.c_str());
					name = PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_PZ).EntityName();
                    break;
                }
			case NNP_NXYZ            :
                {
					name = PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_NX).EntityName();
                    fprintf(fp, "property %s %s\n", type.c_str(), name.c_str());
					name = PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_NY).EntityName();
                    fprintf(fp, "property %s %s\n", type.c_str(), name.c_str());
					name = PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_NZ).EntityName();
                    break;
                }
            case NNP_CRGB            :
                {
					name = PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_CR).EntityName();
                    fprintf(fp, "property %s %s\n", type.c_str(), name.c_str());
					name = PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_CG).EntityName();
                    fprintf(fp, "property %s %s\n", type.c_str(), name.c_str());
					name = PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_CB).EntityName();
                    break;
                }
			case NNP_CRGBA:
				{
					name = PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_CR).EntityName();
					fprintf(fp, "property %s %s\n", type.c_str(), name.c_str());
					name = PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_CG).EntityName();
					fprintf(fp, "property %s %s\n", type.c_str(), name.c_str());
					name = PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_CB).EntityName();
					fprintf(fp, "property %s %s\n", type.c_str(), name.c_str());
					name = PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_CA).EntityName();;
					break;
				}
        }
        fprintf(fp, "property %s %s\n", type.c_str(), name.c_str());
        return true;
    }



	/** PLY Element.
	 *  Define a PLY Element as a collection of properties.
	 */
	class PlyElement
    {
    public:
        
		std::string name;					/**< Name of the elment in the PLY header (for example "vertex" or "face") */
        int cnt;							/**< Number of instances of the elment in the PLY file */
        std::vector<PlyProperty> propVec;   /**< Collection of properties that define the element */

		/** 
		 * Default Constructor
		 */
        PlyElement();

		/** 
		 * Constructor that sets the name, the properties and the number of instances of the element.
		 *
		 * @param *name Name of the element.
		 * @param &prop Vector of properties.
		 * @param nElem Number of instances.
		 */
		PlyElement(const char *name, std::vector<PlyProperty> &prop, int nElem);


		/** 
		 * Parse the input line and add the properties to the element. 
		 * It assumes that the passed line has the folowing structure: "property PLYTYPE PLYELEMENT"
		 * 
		 * @param *line Input line.
		 * @return		If successful returns true. Otherwise, it returns false.
		 */
        bool AddProperty(const char *line);

		 
		/** 
		 * Initialize an element from an header file. 
		 * The first line of the buf is passed into line and must start with the 'element' keyword
		 * 
		 * @param *fp	Pointer to the opened file.
		 * @param *line First buffer line.
		 * @return		If successful returns true. Otherwise, it returns false.
		 */
        bool InitFromHeader(FILE *fp, char *line);
        
		
		/** 
		 * Write the element descriport in the header file. 
		 * 
		 * @param *fp	Pointer to the file.
		 * @return		If successful returns true. Otherwise, it returns false.
		 */
		bool WriteHeader(FILE *fp);
        

		/** 
		 * Skip the element in an Ascii file. 
		 * 
		 * @param *fp	Pointer to the opened file.
		 * @return		If successful returns true. Otherwise, it returns false.
		 */
		bool SkipAsciiElementsInFile(FILE *fp);
        

		/** 
		 * Skip the element in a binary file. 
		 * 
		 * @param *fp	Pointer to the opened file.
		 * @return		If successful returns true. Otherwise, it returns false.
		 */
		bool SkipBinaryElementsInFile(FILE *fp);


		/**
		* Check if the input entity is in the property of the element.
		*
		* @param entity Input entity.
		* @return		If successful returns true. Otherwise, it returns false.
		*/
		bool Contains(NNP_ENTITY entity);
    };


    PlyElement::PlyElement(){}

    PlyElement::PlyElement(const char *_name, std::vector<PlyProperty> &prop, int nElem):name(_name), cnt(nElem), propVec(prop)
    {

    }

    bool PlyElement::InitFromHeader(FILE *fp, char *line)
    {
        char buf[4096];
        for(int i=0;line[i]!=0;++i) buf[i]=tolower(line[i]);
        strtok(buf," \t");    // property
        assert(strstr(buf,"element"));
        char *el = strtok(0," \t\n");
        name = std::string(el);

        sscanf(line,"%*s %*s %i",&cnt);
		//printf("Adding Element '%s' (%i) \n",name.c_str(),cnt);
        fgets(line,4096,fp);
        while(strstr(line,"property"))
        {
            AddProperty(line);
            fgets(line,4096,fp);
        }
        unsigned int mask = 0;
        for (int i = 0; i < propVec.size(); i++)
            mask |= propVec[i].elem;
        std::vector<PlyProperty> compactPropVec;
        for (int i = 0; i < propVec.size(); i++)
        {
            switch (propVec[i].elem)
            {
                case NNP_NX:
                case NNP_NY:
					if ((mask & NNP_NXYZ) != NNP_NXYZ)
                        compactPropVec.push_back(propVec[i]);
                    break;
                case NNP_NZ:
					if ((mask & NNP_NXYZ) != NNP_NXYZ)
                        compactPropVec.push_back(propVec[i]);
                    else
                        compactPropVec.push_back(PlyProperty(propVec[i].type, NNP_NXYZ));
                    break;
                case NNP_PX:
                case NNP_PY:
					if ((mask & NNP_PXYZ) != NNP_PXYZ)
                        compactPropVec.push_back(propVec[i]);
                    break;
                case NNP_PZ:
					if ((mask & NNP_PXYZ) != NNP_PXYZ)
                        compactPropVec.push_back(propVec[i]);
                    else
                        compactPropVec.push_back(PlyProperty(propVec[i].type, NNP_PXYZ));
                    break;
                case NNP_CR:
                case NNP_CG:
					if ((mask & NNP_CRGB) != NNP_CRGB & (mask & NNP_CRGBA) != NNP_CRGBA)
                        compactPropVec.push_back(propVec[i]);
                    break;
                case NNP_CB:
					if ((mask & NNP_CRGB) != NNP_CRGB & (mask & NNP_CRGBA) != NNP_CRGBA)
                        compactPropVec.push_back(propVec[i]);
					else if ((mask & NNP_CRGB) == NNP_CRGB & (mask & NNP_CRGBA) != NNP_CRGBA)
                        compactPropVec.push_back(PlyProperty(propVec[i].type, NNP_CRGB));
                    break;
                case NNP_CA:
					if ((mask & NNP_CRGB) != NNP_CRGB & (mask & NNP_CRGBA) != NNP_CRGBA)
                        compactPropVec.push_back(propVec[i]);
                    else if ((mask & NNP_CRGBA) == NNP_CRGBA)
                        compactPropVec.push_back(PlyProperty(propVec[i].type, NNP_CRGBA));
                    break;
                default:
                    compactPropVec.push_back(propVec[i]);
                    break;
            }
        }

        propVec.clear();
        propVec = compactPropVec;
        return true;
    }


    bool PlyElement::WriteHeader(FILE *fp)
    {
        fprintf(fp, "element %s %d\n", name.c_str(), cnt);
        for(int i = 0; i < propVec.size(); i++)
            propVec[i].WriteHeader(fp);
        return true;
    }


    bool PlyElement::SkipAsciiElementsInFile(FILE *fp)
    {
        char line[4096];
        for(int i = 0; i < this->cnt; ++i)
            fgets(line, 4096, fp);
        return true;
    }


    bool PlyElement::SkipBinaryElementsInFile(FILE *fp)
    {
        for(int i = 0; i < this->cnt; ++i)
            for(int j = 0; j < this->propVec.size(); ++j)
                this->propVec[j].SkipBinaryPropertyInFile(fp);
        return true;
    }


    bool PlyElement::AddProperty(const char *line)
    {
        char buf[128];
        for(int i=0;line[i]!=0;++i) buf[i]=tolower(line[i]);
        strtok(buf," \t");    // property
        char *ty = strtok(0," \t\n"); // float

        PlyType type = NNP_UNKNOWN_TYPE;
        if(strcmp(ty,"float") == 0   || strcmp(ty,"float32") == 0) type = NNP_FLOAT32;
        if(strcmp(ty,"double") == 0 || strcmp(ty,"float64") == 0) type = NNP_FLOAT64;
        if(strcmp(ty,"char") == 0  || strcmp(ty,"int8") == 0)    type = NNP_INT8;
        if(strcmp(ty,"short") == 0  || strcmp(ty,"int16") == 0)   type = NNP_INT16;
        if(strcmp(ty,"int") == 0   || strcmp(ty,"int32") == 0)   type = NNP_INT32;
        if(strcmp(ty,"uchar") == 0 || strcmp(ty,"uint8") == 0 )    type = NNP_UINT8;
        if(strcmp(ty,"ushort") == 0|| strcmp(ty,"uint16") == 0 )   type = NNP_UINT16;
        if(strcmp(ty,"uint") == 0  || strcmp(ty,"uint32") == 0 )   type = NNP_UINT32;
        if(strcmp(ty,"list") == 0 )  {
            char *ty1 = strtok(0," \t\n");
            char *ty2 = strtok(0," \t\n");
            if( (strcmp(ty1,"uchar") == 0 ||strcmp(ty1,"uint8") == 0 ) && (strcmp(ty2,"uint") == 0 || strcmp(ty2,"uint32") == 0) )
				type = NNP_LIST_UINT8_UINT32;
            if( (strcmp(ty1,"char") == 0 || strcmp(ty1,"int8") == 0) && (strcmp(ty2,"uint") == 0|| strcmp(ty2,"uint32") == 0) ) 
				type = NNP_LIST_INT8_UINT32;
            if( (strcmp(ty1,"uchar") == 0 ||strcmp(ty1,"uint8") == 0) && (strcmp(ty2,"int") == 0 || strcmp(ty2,"int32") == 0) ) 
				type = NNP_LIST_UINT8_INT32;
            if( (strcmp(ty1,"char") == 0 || strcmp(ty1,"int8") == 0) && (strcmp(ty2,"int") == 0 || strcmp(ty2,"int32") ) )
				type = NNP_LIST_INT8_INT32;
        }

        assert(type != NNP_UNKNOWN_TYPE);
        char *el = strtok(0," \t\n"); // x

        PlyEntity ent = NNP_UNKNOWN_ENTITY;
        if(strstr(el,"x")) ent = NNP_PX;
        if(strstr(el,"y")) ent = NNP_PY;
        if(strstr(el,"z")) ent = NNP_PZ;
        if(strstr(el,"nx")) ent = NNP_NX;
        if(strstr(el,"ny")) ent = NNP_NY;
        if(strstr(el,"nz")) ent = NNP_NZ;
		if (strstr(el,"red")) ent = NNP_CR;
		if (strstr(el,"diffuse_red")) ent = NNP_CR;
		if (strstr(el,"green")) ent = NNP_CG;
		if (strstr(el,"diffuse_green")) ent = NNP_CG;
		if (strstr(el,"blue")) ent = NNP_CB;
		if (strstr(el,"diffuse_blue")) ent = NNP_CB;
        if(strstr(el,"alpha")) ent = NNP_CA;
		if (strstr(el,"scale")) ent = NNP_SCALE;
		if (strstr(el,"value")) ent = NNP_SCALE;
        if(strstr(el,"density")) ent = NNP_DENSITY;
        if(strstr(el,"radius")) ent = NNP_DENSITY;
		if (strstr(el,"quality")) ent = NNP_QUALITY;
		if (strstr(el,"confidence")) ent = NNP_QUALITY;
        if(strstr(el,"reflectance")) ent = NNP_REFLECTANCE;
        if(strstr(el,"vertex_index")) ent = NNP_VERTEX_LIST;
        if(strstr(el,"vertex_indices")) ent = NNP_VERTEX_LIST;
        //assert(ent != NNP_UNKNOWN_ENTITY);

        propVec.push_back(PlyProperty(type,ent));
        //printf("Adding Property %s %s\n",propVec.back().TypeStr(),propVec.back().EntityStr());

        return true;
    }


	bool PlyElement::Contains(PlyEntity entity)
	{
		for (int i = 0; i < propVec.size(); i++)
		{
			if (propVec[i].elem == entity)
				return true;
		}
		return false;
	}

    
	/** PLY header info.
	 *  Define the data of the PLY header
	 */
	class Info
    {

    public:
        ErrorCode errInfo;					/**< Error code returned by the reading of a PLY file */
        bool binary;						/**< Boolean about the file format (Binary = true, Ascci = false) */
        std::vector<PlyElement> elemVec;	/**< Elements defided in the header */
		bool bigEndian;						/**< Endianess of the binary file */


		/** 
		 * Default Constructor
		 */
        Info();


		/** 
		 * Constructor that reads the header info from a file.
		 *
		 * @param *filename Path of the file to read.
		 */
		Info(const char *filename);


        /** 
		 * Constructor that creates the header info for a mesh.
		 *
		 * @param &vertex	Vertex element object.
		 * @param &face		Face element object.
		 * @param binary	File format (binary = true, ascci = false).
		 */
        Info(PlyElement &vertex, PlyElement &face, bool binary);
		

		/** 
		 * Constructor that creates the header info for a point cloud.
		 *
		 * @param &vertex	Vertex element object.
		 * @param binary	File format (binary = true, ascci = false).
		 */
		Info(PlyElement &vertex, bool binary);


		/** 
		 * Clear the object.
		 */
        void Clear() { errInfo=NNP_OK; }


		/** 
		 * Return the number of vertex instances
		 *
		 * @return The number of vertex instances
		 */
        int GetVertexCount() const;
        

		/** 
		 * Return the number of face instances
		 *
		 * @return The number of face instances
		 */
		int GetFaceCount() const;


		/**
		* Return a reference to the vertex element
		*
		* @return The reference to the vertex element
		*/
		PlyElement* GetVertexElement();


		/**
		* Return a reference to the face element
		*
		* @return The reference to the face element
		*/
		PlyElement* GetFaceElement();
    };


    Info::Info(){}

	Info::Info(const char *filename)
	{
		this->errInfo = NNP_OK;
		FILE *fp=fopen(filename,"r");
      
        if(!fp)
        {
            this->errInfo = NNP_UNABLE_TO_OPEN;
            return;
        }
        char buf[4096];
        fgets(buf,4096,fp);
        if( (strncmp(buf,"PLY",3)!=0) && (strncmp(buf,"ply",3)!=0) )
        {
            this->errInfo = NNP_MISSING_HEADER;
            return;
        }
        fgets(buf,4096,fp);
        if(strncmp(buf,"format",strlen("format"))!=0)
        {
           this->errInfo = NNP_MISSING_FORMAT;
           return;
        }

        if(strstr(buf,"ascii") || strstr(buf,"ASCII"))
        {
            this->binary=false;
        }
        else if(strstr(buf,"binary") || strstr(buf,"BINARY"))
        {
            this->binary=true;
			if (strstr(buf, "binary_big") || strstr(buf, "BINARY_BIG"))
				this->bigEndian = true;
			else
				this->bigEndian = false;
        }
        else
        {
            this->errInfo = NNP_MISSING_FORMAT;
            return;
        }

        fgets(buf,4096,fp);
        while(strncmp(buf,"end_header",strlen("end_header")))
        {
            if(strstr(buf,"comment") || strstr(buf,"COMMENT") )
            {
                fgets(buf,4096,fp);
                continue;
            }
            if(strstr(buf,"element"))
            {
                PlyElement pe;
                pe.InitFromHeader(fp,buf);
                this->elemVec.push_back(pe);
            }
        }
		fclose(fp);
	}

    Info::Info(PlyElement &vertex, PlyElement &face, bool binary)
    {
        elemVec.push_back(vertex);
        elemVec.push_back(face);
        this->binary = binary;
    }

	Info::Info(PlyElement &vertex, bool binary)
    {
        elemVec.push_back(vertex);
        this->binary = binary;
    }

    int Info::GetVertexCount() const
    {
        for (int i = 0; i < elemVec.size(); i++)
        {
            if (elemVec[i].name.compare(std::string("vertex")) == 0)
                return elemVec[i].cnt;
        }
        return -1;
    }

    int Info::GetFaceCount() const
    {
        for (int i = 0; i < elemVec.size(); i++)
        {
            if (elemVec[i].name.compare(std::string("face")) == 0)
                return elemVec[i].cnt;
        }
        return -1;
    }


	PlyElement* Info::GetVertexElement()
	{
		for (int i = 0; i < elemVec.size(); i++)
		{
			if (elemVec[i].name.compare(std::string("vertex")) == 0)
				return &elemVec[i];
		}
		return NULL;
	}


	PlyElement* Info::GetFaceElement()
	{
		for (int i = 0; i < elemVec.size(); i++)
		{
			if (elemVec[i].name.compare(std::string("face")) == 0)
				return &elemVec[i];
		}
		return NULL;
	}


	/**
	 * @cond HIDDEN_SYMBOLS
	 */
	void adjustEndianess(unsigned char* buffer, int typeSize, int count)
	{
		for (int i = 0; i < count; i++)
		{
			int offset = i*typeSize;
			for (int j = 0; j < typeSize/2; j++)
			{
				unsigned char temp = buffer[offset+j];
				buffer[offset+j] = buffer[offset + typeSize - 1 - j];
				buffer[offset + typeSize - 1 - j] = temp;
			}
		}
	}
	/**
	 * @endcond
	 */


	/** Memory descriptor of a vector of properties.
	 *	The class defines how a vector of PlyProperty is saved in memory.
	 *
	 *  @tparam CointainerType	Type of the container of the property 
	 *  @tparam VectorSize		Number of values stored in the property. 
	 *  @tparam ScalarType		Type of the values stored in the property. 
	 */
	template<class CointainerType, int VectorSize, typename ScalarType>
    class VectorDescriptor
    {
        int64_t curPos;			/**< Position of the next property to read or to write. */
        PlyEntity elem;		/**< Ply entity managed by the descriptor. */
        void *base;			/**< Pointer to the memory location that contains the data of the property. */

    public:

		/**
		* Void constructor.
		*
		*/
		VectorDescriptor();


		/** 
		 * Constructor of the descriptor.
		 *
		 * @param _e	Ply entity managed by the descriptor.
		 * @param *_b	Pointer to the memory location that contains the data of the property.
		 */
         VectorDescriptor(PlyEntity _e, void *_b);


		 /** 
		  * Restart the descriptor.
		  */
         void Restart();
         
		 
		 /** 
		  * Read the property data from the binary file if the entity of the property prop is equals to the entity of the descriptor. 
		  *
		  * @param *fp			Pointer to the opened file.
		  * @param &prop		Next PLY property to read from the file.
		  * @param bigEndian	Endianess of the binary data (true = big endian, false = little endian).
		  * @return				If successful returns true. Otherwise, it returns false.
		  */
		 bool ReadElemBinary(FILE *fp, PlyProperty &prop, bool bigEndian);


		 /** 
		  * Read the property data from the ascii file if the entity of the property prop is equals to the entity of the descriptor. 
		  *
		  * @param *fp		Pointer to the opened file.
		  * @param &prop	Next PLY property to read from the file.
		  * @return			If successful returns true. Otherwise, it returns false.
		  */
         bool ReadElemAscii(FILE *fp, PlyProperty &prop);


		 /** 
		  * Write the property data in the binary file if the entity of the property prop is equals to the entity of the descriptor. 
		  *
		  * @param *fp		Pointer to the file.
		  * @param &prop	Next PLY property to write from the file.
		  * @return			If successful returns true. Otherwise, it returns false.
		  */
         bool WriteElemBinary(FILE *fp, PlyProperty &prop);
         
		 
		 /** 
		  * Write the property data in the ascii file if the entity of the property prop is equals to the entity of the descriptor. 
		  *
		  * @param *fp		Pointer to the file.
		  * @param &prop	Next PLY property to write from the file.
		  * @return			If successful returns true. Otherwise, it returns false.
		  */
		 bool WriteElemAscii(FILE *fp, PlyProperty &prop);

    private:

        template<typename C>
        void ReadBinary(FILE *fp, PlyProperty &prop, bool bigEndian);

        template<typename C>
        void ReadAscii(FILE *fp, PlyProperty &prop, const char *format);

        template<typename C>
        void WriteBinary(FILE *fp, PlyProperty &prop);

        template<typename C>
        void WriteAscii(FILE *fp, PlyProperty &prop, const char *format);

    };


	template<class CointainerType, int VectorSize, typename ScalarType>
	VectorDescriptor<CointainerType, VectorSize, ScalarType>::VectorDescriptor()
	{

	}

    template<class CointainerType, int VectorSize, typename ScalarType>
    VectorDescriptor<CointainerType, VectorSize, ScalarType>::VectorDescriptor(PlyEntity _e, void *_b):curPos(0),elem(_e),base(_b)
    {

    }

    template<class CointainerType, int VectorSize, typename ScalarType>
    void VectorDescriptor<CointainerType, VectorSize, ScalarType>::Restart()
    {
        this->curPos=0;
    }

    template<class ContainerType, int VectorSize, typename ScalarType>
    template<typename C>
    void VectorDescriptor<ContainerType, VectorSize, ScalarType>::ReadBinary(FILE *fp, PlyProperty &prop, bool bigEndian)
    {
        int count = 1;
        if (prop.elem == NNP_CRGB || prop.elem == NNP_NXYZ || prop.elem == NNP_PXYZ)
            count = 3;
        else if (prop.elem == NNP_CRGBA)
            count = 4;
        int typeSize = prop.TypeSize();
        unsigned char* buffer = new unsigned char[count*typeSize];
        fread(buffer, typeSize, count, fp);
		if (typeSize > 1 && bigEndian)
			adjustEndianess(buffer, typeSize, count);
        if (prop.type == NNP_LIST_UINT8_UINT32 || prop.type == NNP_LIST_UINT8_INT32)
        {
            count = buffer[0];
            delete[] buffer;
            buffer = new unsigned char[count*4];
            fread(buffer, 4, count, fp);
			if (bigEndian)
				adjustEndianess(buffer, 4, count);
        }
        else if (prop.type == NNP_LIST_INT8_UINT32 || prop.type == NNP_LIST_INT8_INT32)
        {
            count = char(buffer[0]);
            delete[] buffer;
            buffer = new unsigned char[count*4];
            fread(buffer, 4, count, fp);
			if (bigEndian)
				adjustEndianess(buffer, 4, count);
        }
        C* temp = (C*)buffer;
        
		unsigned char* baseProp = (unsigned char*)base + this->curPos*sizeof(ContainerType);
		for (int i = 0; i < std::min(VectorSize, count); i++)
			*(ScalarType *)(baseProp + i*sizeof(ScalarType)) = ScalarType(temp[i]);
		
		++(this->curPos);
        delete[] buffer;
    }

	template<class ContainerType, int VectorSize, typename ScalarType>
    bool VectorDescriptor<ContainerType, VectorSize, ScalarType>::ReadElemBinary(FILE *fp, PlyProperty &prop, bool bigEndian)
    {
        if (prop.elem != elem)
            return false;
        switch(prop.type)
        {
            case NNP_INT8:				this->ReadBinary<char>(fp, prop, bigEndian); break;
            case NNP_UINT8:				this->ReadBinary<unsigned char>(fp, prop, bigEndian); break;
            case NNP_INT16:				this->ReadBinary<short>(fp, prop, bigEndian); break;
            case NNP_UINT16:			this->ReadBinary<unsigned short>(fp, prop, bigEndian); break;
            case NNP_FLOAT32:			this->ReadBinary<float>(fp, prop, bigEndian); break;
            case NNP_INT32:				this->ReadBinary<int>(fp, prop, bigEndian); break;
            case NNP_UINT32:			this->ReadBinary<unsigned int>(fp, prop, bigEndian); break;
            case NNP_FLOAT64:			this->ReadBinary<double>(fp, prop, bigEndian); break;
            case NNP_LIST_UINT8_UINT32:
            case NNP_LIST_INT8_UINT32:	this->ReadBinary<unsigned int>(fp, prop, bigEndian); break;
            case NNP_LIST_UINT8_INT32:
            case NNP_LIST_INT8_INT32:	this->ReadBinary<int>(fp, prop, bigEndian); break;
        }
        return true;
    }



    template<class ContainerType, int VectorSize, typename ScalarType>
    template<typename C>
    void VectorDescriptor<ContainerType, VectorSize, ScalarType>::ReadAscii(FILE *fp, PlyProperty &prop, const char *format)
    {
        int count = 1;
        if (prop.elem == NNP_CRGB || prop.elem == NNP_NXYZ || prop.elem == NNP_PXYZ)
            count = 3;
        else if (prop.elem == NNP_CRGBA)
            count = 4;

        if (prop.type == NNP_LIST_UINT8_UINT32 || prop.type == NNP_LIST_UINT8_INT32)
        {
            unsigned char listSize;
            fscanf(fp, "%u", &listSize);
            count = listSize;
        }
        else if (prop.type == NNP_LIST_INT8_UINT32 || prop.type == NNP_LIST_INT8_INT32)
        {
            char listSize;
            fscanf(fp, "%d", &listSize);
            count = listSize;
        }

        C* temp = new C[count];
        for (int i = 0; i < count ; i++)
            fscanf(fp, format, &temp[i]);
       
		unsigned char* baseProp = (unsigned char*)base + this->curPos*sizeof(ContainerType);
		for (int i = 0; i < std::min(VectorSize,count); i++)
			*(ScalarType *)(baseProp + i*sizeof(ScalarType)) = ScalarType(temp[i]);	
		
		delete[] temp;
        ++(this->curPos);
    }



    template<class ContainerType, int VectorSize, typename ScalarType>
    bool VectorDescriptor<ContainerType, VectorSize, ScalarType>::ReadElemAscii(FILE *fp, PlyProperty &prop)
    {
        if (prop.elem != elem)
            return false;
        switch(prop.type)
        {
            case NNP_INT8:				this->ReadAscii<int>(fp, prop, "%d"); break;
            case NNP_UINT8:				this->ReadAscii<unsigned int>(fp, prop, "%u"); break;
            case NNP_INT16:				this->ReadAscii<short>(fp, prop, "%d"); break;
            case NNP_UINT16:			this->ReadAscii<unsigned short>(fp, prop, "%u"); break;
            case NNP_FLOAT32:			this->ReadAscii<float>(fp, prop, "%f"); break;
            case NNP_INT32:				this->ReadAscii<int>(fp, prop, "%d"); break;
            case NNP_UINT32:			this->ReadAscii<unsigned int>(fp, prop, "%u"); break;
            case NNP_FLOAT64:			this->ReadAscii<double>(fp, prop, "%f"); break;
            case NNP_LIST_UINT8_UINT32:
            case NNP_LIST_INT8_UINT32:	this->ReadAscii<unsigned int>(fp, prop, "%u"); break;
            case NNP_LIST_UINT8_INT32:
            case NNP_LIST_INT8_INT32:	this->ReadAscii<int>(fp, prop, "%d"); break;
        }
        return true;
    }



    template<class ContainerType, int VectorSize, typename ScalarType>
    template<typename C>
    void VectorDescriptor<ContainerType, VectorSize, ScalarType>::WriteBinary(FILE *fp, PlyProperty &prop)
    {
        int count = 1;
        if (prop.elem == NNP_CRGB || prop.elem == NNP_NXYZ || prop.elem == NNP_PXYZ)
            count = 3;
        else if (prop.elem == NNP_CRGBA)
            count = 4;

		C data[VectorSize];
		
		if (prop.type == NNP_LIST_UINT8_UINT32 || prop.type == NNP_LIST_UINT8_INT32)
        {
            unsigned char listSize = (unsigned char) VectorSize;
            fwrite(&listSize, sizeof(unsigned char), 1, fp);
            count = VectorSize;
        }
        else if (prop.type == NNP_LIST_INT8_UINT32 || prop.type == NNP_LIST_INT8_INT32)
        {
            char listSize = (char) VectorSize;
            fwrite(&listSize, sizeof(char), 1, fp);
            count = VectorSize;
        }
		C temp = 0;
		unsigned char* baseProp = (unsigned char*)base + this->curPos*sizeof(ContainerType);
		for (int i = 0; i < std::min(VectorSize, count); i++)
			data[i] = (C)(*(ScalarType*)(baseProp + i*sizeof(ScalarType)));
		fwrite(data, sizeof(C), std::min(VectorSize, count), fp);
		for (int i = 0; i < (count - VectorSize); i++)
			fwrite(&temp, sizeof(C), 1, fp);        
        ++(this->curPos);
    }


    template<class ContainerType, int VectorSize, typename ScalarType>
    bool VectorDescriptor<ContainerType, VectorSize, ScalarType>::WriteElemBinary(FILE *fp, PlyProperty &prop)
    {
        if (prop.elem != elem)
            return false;
        switch(prop.type)
        {
            case NNP_INT8:				this->WriteBinary<char>(fp, prop); break;
            case NNP_UINT8:				this->WriteBinary<unsigned char>(fp, prop); break;
            case NNP_INT16:				this->WriteBinary<short>(fp, prop); break;
            case NNP_UINT16:			this->WriteBinary<unsigned short>(fp, prop); break;
            case NNP_FLOAT32:			this->WriteBinary<float>(fp, prop); break;
            case NNP_INT32:				this->WriteBinary<int>(fp, prop); break;
            case NNP_UINT32:			this->WriteBinary<unsigned int>(fp, prop); break;
            case NNP_FLOAT64:			this->WriteBinary<double>(fp, prop); break;
            case NNP_LIST_UINT8_UINT32:
            case NNP_LIST_INT8_UINT32:	this->WriteBinary<unsigned int>(fp, prop); break;
            case NNP_LIST_UINT8_INT32:
            case NNP_LIST_INT8_INT32:	this->WriteBinary<int>(fp, prop); break;
        }
        return true;
    }


    template<class ContainerType, int VectorSize, typename ScalarType>
    template<typename C>
    void VectorDescriptor<ContainerType, VectorSize, ScalarType>::WriteAscii(FILE *fp, PlyProperty &prop, const char *format)
    {
        int count = 1;
        if (prop.elem == NNP_CRGB || prop.elem == NNP_NXYZ || prop.elem == NNP_PXYZ)
            count = 3;
        else if (prop.elem == NNP_CRGBA)
            count = 4;

        if (prop.type == NNP_LIST_UINT8_UINT32 || prop.type == NNP_LIST_UINT8_INT32)
        {
            fprintf(fp, "%u ", (unsigned char)(VectorSize));
            count = VectorSize;
        }
        else if (prop.type == NNP_LIST_INT8_UINT32 || prop.type == NNP_LIST_INT8_INT32)
        {
            fprintf(fp, "%d ", (char)(VectorSize));
            count = VectorSize;
        }

		unsigned char* baseProp = (unsigned char*)base + this->curPos*sizeof(ContainerType);
        for (int i = 0; i < std::min(VectorSize, count); i++)
            fprintf(fp, format, (C)(*(ScalarType*)(baseProp + i*sizeof(ScalarType))));
        for (int i = 0; i < (count - VectorSize); i++)
            fprintf(fp, format, (C)(0));
        ++(this->curPos);
    }

    template<class ContainerType, int VectorSize, typename ScalarType>
    bool VectorDescriptor<ContainerType, VectorSize, ScalarType>::WriteElemAscii(FILE *fp, PlyProperty& prop)
    {
        if (prop.elem != elem)
            return false;

        switch(prop.type)
        {
            case NNP_INT8:				this->WriteAscii<char>(fp, prop, "%d "); break;
            case NNP_UINT8:				this->WriteAscii<unsigned char>(fp, prop, "%u "); break;
            case NNP_INT16:				this->WriteAscii<short>(fp, prop, "%d "); break;
            case NNP_UINT16:			this->WriteAscii<unsigned short>(fp, prop, "%u "); break;
            case NNP_FLOAT32:			this->WriteAscii<float>(fp, prop, "%f "); break;
            case NNP_INT32:				this->WriteAscii<int>(fp, prop, "%d "); break;
            case NNP_UINT32:			this->WriteAscii<unsigned int>(fp, prop, "%u "); break;
            case NNP_FLOAT64:			this->WriteAscii<double>(fp, prop, "%f "); break;
            case NNP_LIST_UINT8_UINT32:
            case NNP_LIST_INT8_UINT32:	this->WriteAscii<unsigned int>(fp, prop, "%u "); break;
            case NNP_LIST_UINT8_INT32:
            case NNP_LIST_INT8_INT32:	this->WriteAscii<int>(fp, prop, "%d "); break;
        }
        return true;
    }




	/**
	 * @cond HIDDEN_SYMBOLS
	 */
    template <class TupleType>
    void ReadBinaryElement(TupleType &adaptor, PlyElement &elem, FILE *fp, bool bigEndian)
    {
        for (int i = 0 ; i < elem.cnt; i++)
        {
            for (int j = 0; j < elem.propVec.size(); j++)
            {
                PlyProperty& prop = elem.propVec[j];
               	if (!TupleForEach(adaptor, fp, prop, bigEndian, SizeT<0>()))
                    prop.SkipBinaryPropertyInFile(fp);
            }
        }
    }


    template <class TupleType>
    void ReadAsciiElement(TupleType &adaptor, PlyElement &elem, FILE *fp)
    {
        for (int i = 0 ; i < elem.cnt; i++)
        {
            for (int j = 0; j < elem.propVec.size(); j++)
            {
                PlyProperty& prop = elem.propVec[j];
                if (!TupleForEach(adaptor, fp, prop, false, SizeT<1>()))
                    prop.SkipAsciiPropertyInFile(fp);
            }
        }
    }
	/**
	 * @endcond
	 */



    /**
	 * Read a point cloud from a PLY file.
	 *
	 * @tparam VertexAdaptorTuple	Type that defines the management of the vertex data in memory
	 * 
	 * @param filename				Path to the file to read
	 * @param vertexAdaptor			std::tuple that defines how to manage the vertex data in memory
	 */
    template <class VertexAdaptorTuple>
    bool OpenPointCloud(const char *filename, VertexAdaptorTuple vertexAdaptor)
    {
		nanoply::Info info(filename);
		if (info.errInfo != NNP_OK)
			return false;
        FILE *fp=fopen(filename,"rb");
        if(!fp)
        {
            return false;
        }
        char buf[4096];
        do
        {
            fgets(buf,4096,fp);
        }
        while(strncmp(buf,"end_header",strlen("end_header")) );
        // Now start the real reading!

        if (info.binary)
        {
            for(int i = 0; i < info.elemVec.size();++i)
                if (strcmp(info.elemVec[i].name.c_str(), "vertex") == 0)
					ReadBinaryElement<VertexAdaptorTuple>(vertexAdaptor, info.elemVec[i], fp, info.bigEndian);
                else
                    info.elemVec[i].SkipBinaryElementsInFile(fp);
        }
        else
        {
            for(int i = 0; i < info.elemVec.size();++i)
                if (strcmp(info.elemVec[i].name.c_str(), "vertex") == 0)
                    ReadAsciiElement<VertexAdaptorTuple>(vertexAdaptor, info.elemVec[i], fp);
                else
                    info.elemVec[i].SkipAsciiElementsInFile(fp);
        }
        return true;
    }


    /**
	 * Read a mesh from a PLY file.
	 *
	 * @tparam VertexAdaptorTuple	Type that defines the management of the vertex data in memory
	 * @tparam FaceAdaptorTuple		Type that defines the management of the face data in memory
	 * 
	 * @param filename				Path to the file to read
	 * @param vertexAdaptor			std::tuple that defines how to manage the vertex data in memory
	 * @param faceAdaptor			std::tuple that defines how to manage the face data in memory
	 */
    template <class VertexAdaptorTuple, class FaceAdaptorTuple>
    bool OpenMesh(const char *filename, VertexAdaptorTuple vertexAdaptor, FaceAdaptorTuple faceAdaptor)
    {
		nanoply::Info info(filename);
		if (info.errInfo != NNP_OK)
			return false;
        FILE *fp=fopen(filename,"rb");
        if(!fp)
        {
            return false;
        }
        char buf[4096];
        do
        {
            fgets(buf,4096,fp);
        }
        while(strncmp(buf,"end_header",strlen("end_header")) );
        // Now start the real reading!

        if (info.binary)
        {
            for(int i = 0; i < info.elemVec.size();++i)
            {
                if (strcmp(info.elemVec[i].name.c_str(), "vertex") == 0)
					ReadBinaryElement<VertexAdaptorTuple>(vertexAdaptor, info.elemVec[i], fp, info.bigEndian);
                else if (strcmp(info.elemVec[i].name.c_str(), "face") == 0)
					ReadBinaryElement<FaceAdaptorTuple>(faceAdaptor, info.elemVec[i], fp, info.bigEndian);
                else
                    info.elemVec[i].SkipBinaryElementsInFile(fp);
            }
        }
        else
        {
            for(int i = 0; i < info.elemVec.size();++i)
            {
                if (strcmp(info.elemVec[i].name.c_str(), "vertex") == 0)
                    ReadAsciiElement<VertexAdaptorTuple>(vertexAdaptor, info.elemVec[i], fp);
                else if (strcmp(info.elemVec[i].name.c_str(), "face") == 0)
                    ReadAsciiElement<FaceAdaptorTuple>(faceAdaptor, info.elemVec[i], fp);
                else
                    info.elemVec[i].SkipAsciiElementsInFile(fp);
            }
        }
    }



	/**
	 * @cond HIDDEN_SYMBOLS
	 */
    template <class TupleType>
    void WriteBinaryElement(TupleType &adaptor, PlyElement &elem, FILE *fp)
    {
        for (int i = 0 ; i < elem.cnt; i++)
        {
            for (int j = 0; j < elem.propVec.size(); j++)
            {
                PlyProperty& prop = elem.propVec[j];
                TupleForEach(adaptor, fp, prop, false, SizeT<2>());
            }
        }
    }


    template <class TupleType>
    void WriteAsciiElement(TupleType &adaptor, PlyElement &elem, FILE *fp)
    {
        for (int i = 0 ; i < elem.cnt; i++)
        {
            for (int j = 0; j < elem.propVec.size(); j++)
            {
                PlyProperty& prop = elem.propVec[j];
				TupleForEach(adaptor, fp, prop, false, SizeT<3>());
            }
            fprintf(fp, "\n");
        }
    }

	/**
	 * @endcond
	 */


    /**
	 * Save a point cloud from a PLY file.
	 *
	 * @tparam VertexAdaptorTuple	Type that defines the management of the vertex data in memory
	 * 
	 * @param filename				Path to the file to save
	 * @param vertexAdaptor			std::tuple that defines how to manage the vertex data in memory
	 * @param info					Info to saved in the PLY header
	 */
	template <class VertexAdaptorTuple>
    bool SavePointCloud(const char *filename, VertexAdaptorTuple vertexAdaptor, nanoply::Info &info)
    {
        FILE *fp=fopen(filename,"wb");
        if(!fp)
        {
            return false;
        }
        fprintf(fp, "ply\n");
        if (info.binary)
            fprintf(fp, "format binary_little_endian 1.0\n");
        else
            fprintf(fp, "format ascii 1.0\n");

        for (int i = 0; i < info.elemVec.size(); i++)
        {
            if (strcmp(info.elemVec[i].name.c_str(), "vertex") != 0)
            {
                info.elemVec[i].cnt = 0;
                info.elemVec[i].WriteHeader(fp);
            }
            else
                info.elemVec[i].WriteHeader(fp);
        }

        fprintf(fp, "end_header\n");

        if (info.binary)
        {
            for(int i = 0; i < info.elemVec.size();++i)
                if (strcmp(info.elemVec[i].name.c_str(), "vertex") == 0)
                    WriteBinaryElement<VertexAdaptorTuple>(vertexAdaptor, info.elemVec[i], fp);
        }
        else
        {
            for(int i = 0; i < info.elemVec.size();++i)
                if (strcmp(info.elemVec[i].name.c_str(), "vertex") == 0)
                    WriteAsciiElement<VertexAdaptorTuple>(vertexAdaptor, info.elemVec[i], fp);
        }

        fclose(fp);
        return true;
    }


    /**
	 * Save a mesh from a PLY file.
	 *
	 * @tparam VertexAdaptorTuple	Type that defines the management of the vertex data in memory
	 * @tparam FaceAdaptorTuple		Type that defines the management of the face data in memory
	 * 
	 * @param filename				Path to the file to save
	 * @param vertexAdaptor			std::tuple that defines how to manage the vertex data in memory
	 * @param faceAdaptor			std::tuple that defines how to manage the face data in memory
	 * @param info					Info to saved in the PLY header
	 */
    template <class VertexAdaptorTuple, class FaceAdaptorTuple>
    bool SaveMesh(const char *filename, VertexAdaptorTuple vertexAdaptor, FaceAdaptorTuple faceAdaptor, nanoply::Info &info)
    {
        FILE *fp=fopen(filename,"wb");
        if(!fp)
        {
            return false;
        }
        fprintf(fp, "ply\n");
        if (info.binary)
            fprintf(fp, "format binary_little_endian 1.0\n");
        else
            fprintf(fp, "format ascii 1.0\n");

        for (int i = 0; i < info.elemVec.size(); i++)
        {
            if (strcmp(info.elemVec[i].name.c_str(), "vertex") != 0 && strcmp(info.elemVec[i].name.c_str(), "face") != 0)
            {
                info.elemVec[i].cnt = 0;
                info.elemVec[i].WriteHeader(fp);
            }
            else
                info.elemVec[i].WriteHeader(fp);
        }

        fprintf(fp, "end_header\n");

        if (info.binary)
        {
            for(int i = 0; i < info.elemVec.size();++i)
            {
                if (strcmp(info.elemVec[i].name.c_str(), "vertex") == 0)
                    WriteBinaryElement<VertexAdaptorTuple>(vertexAdaptor, info.elemVec[i], fp);
                else if (strcmp(info.elemVec[i].name.c_str(), "face") == 0)
                    WriteBinaryElement<FaceAdaptorTuple>(faceAdaptor, info.elemVec[i], fp);
            }
        }
        else
        {
            for(int i = 0; i < info.elemVec.size();++i)
            {
                if (strcmp(info.elemVec[i].name.c_str(), "vertex") == 0)
                    WriteAsciiElement<VertexAdaptorTuple>(vertexAdaptor, info.elemVec[i], fp);
                else if (strcmp(info.elemVec[i].name.c_str(), "face") == 0)
                    WriteAsciiElement<FaceAdaptorTuple>(faceAdaptor, info.elemVec[i], fp);
            }
        }

        fclose(fp);
        return true;
    }



	/**
	 * @cond HIDDEN_SYMBOLS
	 */
   	template < typename TupleType, size_t ActionType>
    inline bool TupleForEach( TupleType &tuple, FILE *fp, PlyProperty &pro, bool bigEndian, SizeT<ActionType> a)
    {
		return TupleForEach( tuple, pro, fp, bigEndian, SizeT<std::tuple_size<TupleType>::value>(), a);
    }

    template < typename TupleType, size_t ActionType>
    inline bool TupleForEach( TupleType &tuple, PlyProperty &pro, FILE *fp, bool bigEndian, SizeT<0> t, SizeT<ActionType> a) {return false; }

    template < typename TupleType, size_t N, size_t ActionType>
    inline bool TupleForEach( TupleType &tuple, PlyProperty &pro, FILE *fp, bool bigEndian, SizeT<N> t, SizeT<ActionType> a)
    {
        typename std::tuple_element<N - 1, TupleType>::type &dataDescr = std::get<N -1>(tuple);
		if (ActionType == 0)
		{
			if (dataDescr.ReadElemBinary(fp, pro, bigEndian))
				return true;
		}
		else if (ActionType == 1)
		{
			if (dataDescr.ReadElemAscii(fp, pro))
				return true;
		}
		else if (ActionType == 2)
		{
			if (dataDescr.WriteElemBinary(fp, pro))
				return true;
		}
		else if (ActionType == 3)
		{
			if (dataDescr.WriteElemAscii(fp, pro))
				return true;
		}
        return TupleForEach( tuple, pro, fp, bigEndian, SizeT<N-1>(), a);
    }

	/**
	 * @endcond
	 */

} // end namespace nanoply
#endif // NANOPLY_HPP
