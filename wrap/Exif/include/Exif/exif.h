#ifndef EXIF_H
#define EXIF_H

#include <vcg/math/camera.h>
#include <vector>
using namespace std;

template <class T>
class Exif{
public:
  Exif(){
    sections.resize(0);
    motorolaOrder = false;
  }

  bool readJPEG(const char *filename){
    FILE *fp = fopen(filename, "rb"); // Unix ignores 'b', windows needs it.
    //parse the marker stream until SOS or EOI is seen
    if(!fp || !readJpegSections(fp)) return false;
    fclose(fp);
    return true;
  }

  //Compute the CCD width, in millimeters
  // Note: With some cameras, its not possible to compute this correctly because
  // they don't adjust the indicated focal plane resolution units when using less
  // than maximum resolution, so the CCDWidth value comes out too small. Nothing
  // that Jhad can do about it - its a camera problem.
  T CCDwidthMm() const {
    if(!FocalplaneXRes) return -1;
    return T(ExifImageWidth * FocalplaneUnits / FocalplaneXRes);
  }

  // Compute 35 mm equivalent focal length based on sensor geometry if we haven't
  // already got it explicitly from a tag.
  T FocalLength35mmEquiv() const {
    if( intrinsics.FocalMm && !FocalLength35mmEquiv )
      FocalLength35mmEquiv = T(ImageInfo.FocalLength/ImageInfo.CCDWidth*36 + 0.5);
    return FocalLength35mmEquiv;
  }

  vcg::Camera<T> vcgCamera() const {
    vcg::Camera<T> intrinsics;
    return intrinsics;
  }

private:
// JPEG markers consist of one or more 0xFF bytes, followed by a marker
// code byte (which is not an FF). Here are the marker codes of interest
// in this program.
#define M_SOF0  0xC0          // Start Of Frame N
#define M_SOF1  0xC1          // N indicates which compression process
#define M_SOF2  0xC2          // Only SOF0-SOF2 are now in common use
#define M_SOF3  0xC3
#define M_SOF5  0xC5          // NB: codes C4 and CC are NOT SOF markers
#define M_SOF6  0xC6
#define M_SOF7  0xC7
#define M_SOF9  0xC9
#define M_SOF10 0xCA
#define M_SOF11 0xCB
#define M_SOF13 0xCD
#define M_SOF14 0xCE
#define M_SOF15 0xCF
#define M_SOI   0xD8          // Start Of Image (beginning of datastream)
#define M_EOI   0xD9          // End Of Image (end of datastream)
#define M_SOS   0xDA          // Start Of Scan (begins compressed data)
#define M_JFIF  0xE0          // Jfif marker
#define M_EXIF  0xE1          // Exif marker.  Also used for XMP data!
#define M_XMP   0x10E1        // Not a real tag (same value in file as Exif!)
#define M_COM   0xFE          // COMment 
#define M_DQT   0xDB
#define M_DHT   0xC4
#define M_DRI   0xDD
#define M_IPTC  0xED          // IPTC marker

// Exif format descriptor stuff
#define MAX_DATE_COPIES 10
#ifdef _WIN32
    #define PATH_MAX _MAX_PATH
    #define SLASH '\\'
#else
    #define SLASH '/'
#endif
const int *BytesPerFormat;
#define NUM_FORMATS 12
#define FMT_BYTE       1 
#define FMT_STRING     2
#define FMT_USHORT     3
#define FMT_ULONG      4
#define FMT_URATIONAL  5
#define FMT_SBYTE      6
#define FMT_UNDEFINED  7
#define FMT_SSHORT     8
#define FMT_SLONG      9
#define FMT_SRATIONAL 10
#define FMT_SINGLE    11
#define FMT_DOUBLE    12

// Describes tag values
#define TAG_INTEROP_INDEX          0x0001
#define TAG_INTEROP_VERSION        0x0002
#define TAG_IMAGE_WIDTH            0x0100
#define TAG_IMAGE_LENGTH           0x0101
#define TAG_BITS_PER_SAMPLE        0x0102
#define TAG_COMPRESSION            0x0103
#define TAG_PHOTOMETRIC_INTERP     0x0106
#define TAG_FILL_ORDER             0x010A
#define TAG_DOCUMENT_NAME          0x010D
#define TAG_IMAGE_DESCRIPTION      0x010E
#define TAG_MAKE                   0x010F
#define TAG_MODEL                  0x0110
#define TAG_SRIP_OFFSET            0x0111
#define TAG_ORIENTATION            0x0112
#define TAG_SAMPLES_PER_PIXEL      0x0115
#define TAG_ROWS_PER_STRIP         0x0116
#define TAG_STRIP_BYTE_COUNTS      0x0117
#define TAG_X_RESOLUTION           0x011A
#define TAG_Y_RESOLUTION           0x011B
#define TAG_PLANAR_CONFIGURATION   0x011C
#define TAG_RESOLUTION_UNIT        0x0128
#define TAG_TRANSFER_FUNCTION      0x012D
#define TAG_SOFTWARE               0x0131
#define TAG_DATETIME               0x0132
#define TAG_ARTIST                 0x013B
#define TAG_WHITE_POINT            0x013E
#define TAG_PRIMARY_CHROMATICITIES 0x013F
#define TAG_TRANSFER_RANGE         0x0156
#define TAG_JPEG_PROC              0x0200
#define TAG_THUMBNAIL_OFFSET       0x0201
#define TAG_THUMBNAIL_LENGTH       0x0202
#define TAG_Y_CB_CR_COEFFICIENTS   0x0211
#define TAG_Y_CB_CR_SUB_SAMPLING   0x0212
#define TAG_Y_CB_CR_POSITIONING    0x0213
#define TAG_REFERENCE_BLACK_WHITE  0x0214
#define TAG_RELATED_IMAGE_WIDTH    0x1001
#define TAG_RELATED_IMAGE_LENGTH   0x1002
#define TAG_CFA_REPEAT_PATTERN_DIM 0x828D
#define TAG_CFA_PATTERN1           0x828E
#define TAG_BATTERY_LEVEL          0x828F
#define TAG_COPYRIGHT              0x8298
#define TAG_EXPOSURETIME           0x829A
#define TAG_FNUMBER                0x829D
#define TAG_IPTC_NAA               0x83BB
#define TAG_EXIF_OFFSET            0x8769
#define TAG_INTER_COLOR_PROFILE    0x8773
#define TAG_EXPOSURE_PROGRAM       0x8822
#define TAG_SPECTRAL_SENSITIVITY   0x8824
#define TAG_GPSINFO                0x8825
#define TAG_ISO_EQUIVALENT         0x8827
#define TAG_OECF                   0x8828
#define TAG_EXIF_VERSION           0x9000
#define TAG_DATETIME_ORIGINAL      0x9003
#define TAG_DATETIME_DIGITIZED     0x9004
#define TAG_COMPONENTS_CONFIG      0x9101
#define TAG_CPRS_BITS_PER_PIXEL    0x9102
#define TAG_SHUTTERSPEED           0x9201
#define TAG_APERTURE               0x9202
#define TAG_BRIGHTNESS_VALUE       0x9203
#define TAG_EXPOSURE_BIAS          0x9204
#define TAG_MAXAPERTURE            0x9205
#define TAG_SUBJECT_DISTANCE       0x9206
#define TAG_METERING_MODE          0x9207
#define TAG_LIGHT_SOURCE           0x9208
#define TAG_FLASH                  0x9209
#define TAG_FOCALLENGTH            0x920A
#define TAG_MAKER_NOTE             0x927C
#define TAG_USERCOMMENT            0x9286
#define TAG_SUBSEC_TIME            0x9290
#define TAG_SUBSEC_TIME_ORIG       0x9291
#define TAG_SUBSEC_TIME_DIG        0x9292
#define TAG_WINXP_TITLE            0x9c9b // Windows XP - not part of exif standard.
#define TAG_WINXP_COMMENT          0x9c9c // Windows XP - not part of exif standard.
#define TAG_WINXP_AUTHOR           0x9c9d // Windows XP - not part of exif standard.
#define TAG_WINXP_KEYWORDS         0x9c9e // Windows XP - not part of exif standard.
#define TAG_WINXP_SUBJECT          0x9c9f // Windows XP - not part of exif standard.
#define TAG_FLASH_PIX_VERSION      0xA000
#define TAG_COLOR_SPACE            0xA001
#define TAG_EXIF_IMAGEWIDTH        0xA002
#define TAG_EXIF_IMAGELENGTH       0xA003
#define TAG_RELATED_AUDIO_FILE     0xA004
#define TAG_INTEROP_OFFSET         0xA005
#define TAG_FLASH_ENERGY           0xA20B
#define TAG_SPATIAL_FREQ_RESP      0xA20C
#define TAG_FOCAL_PLANE_XRES       0xA20E
#define TAG_FOCAL_PLANE_YRES       0xA20F
#define TAG_FOCAL_PLANE_UNITS      0xA210
#define TAG_SUBJECT_LOCATION       0xA214
#define TAG_EXPOSURE_INDEX         0xA215
#define TAG_SENSING_METHOD         0xA217
#define TAG_FILE_SOURCE            0xA300
#define TAG_SCENE_TYPE             0xA301
#define TAG_CFA_PATTERN            0xA302
#define TAG_CUSTOM_RENDERED        0xA401
#define TAG_EXPOSURE_MODE          0xA402
#define TAG_WHITEBALANCE           0xA403
#define TAG_DIGITALZOOMRATIO       0xA404
#define TAG_FOCALLENGTH_35MM       0xA405
#define TAG_SCENE_CAPTURE_TYPE     0xA406
#define TAG_GAIN_CONTROL           0xA407
#define TAG_CONTRAST               0xA408
#define TAG_SATURATION             0xA409
#define TAG_SHARPNESS              0xA40A
#define TAG_DISTANCE_RANGE         0xA40C

  typedef unsigned int uint;
  typedef unsigned char uchar;
  static const uint MAX_COMMENT_SIZE = 2000;

  typedef struct Section{
    int type;
    uint size;
    uchar *data;
  }Section;
  vector<Section> sections;

  bool motorolaOrder;

  // This structure stores Exif header image elements in a simple manner
  // Used to store camera data as extracted from the various ways that it can be
  // stored in an exif header
  typedef struct {
      char  FileName     [PATH_MAX+1];
      time_t FileDateTime;
      unsigned FileSize;
      char  CameraMake   [32];
      char  CameraModel  [40];
      char  DateTime     [20];
      int   Height, Width;
      int   Orientation;
      int   IsColor;
      int   Process;
      int   FlashUsed;
      float FocalLength;
      float ExposureTime;
      float ApertureFNumber;
      float Distance;
      float CCDWidth;
      float ExposureBias;
      float DigitalZoomRatio;
      int   FocalLength35mmEquiv; // Exif 2.2 tag - usually not present.
      int   Whitebalance;
      int   MeteringMode;
      int   ExposureProgram;
      int   ExposureMode;
      int   ISOequivalent;
      int   LightSource;
      int   DistanceRange;

      char  Comments[MAX_COMMENT_SIZE];
      int   CommentWidchars; // If nonzer, widechar comment, indicates number of chars.

      unsigned ThumbnailOffset;          // Exif offset to thumbnail
      unsigned ThumbnailSize;            // Size of thumbnail.
      unsigned LargestExifOffset;        // Last exif data referenced (to check if thumbnail is at end)

      char  ThumbnailAtEnd;              // Exif header ends with the thumbnail
                                         // (we can only modify the thumbnail if its at the end)
      int   ThumbnailSizeOffset;

      int  DateTimeOffsets[MAX_DATE_COPIES];
      int  numDateTimeTags;

      int GpsInfoPresent;
      char GpsLat[31];
      char GpsLong[31];
      char GpsAlt[20];
  }ImageInfo_t;

  ImageInfo_t ImageInfo;

  uchar *DirWithThumbnailPtrs;
  T FocalplaneXRes;
  T FocalplaneUnits;
  int ExifImageWidth;
  void * OrientationPtr[2];
  int    OrientationNumFormat[2];
  int NumOrientations;

  // Get 16 bits motorola order (always) for jpeg header stuff.
  int Get16m(const void *Short){ return (((uchar*)Short)[0] << 8) | ((uchar*)Short)[1]; }
  // Convert a 16 bit unsigned value from file's native byte order
  int Get16u(void *Short){
    return (motorolaOrder)?
      (((uchar *)Short)[0] << 8) | ((uchar *)Short)[1] :
      (((uchar *)Short)[1] << 8) | ((uchar *)Short)[0] ;
  }
  // Convert a 32 bit unsigned value from file's native byte order
  unsigned Get32u(void *Long){ return (unsigned)Get32s(Long) & 0xffffffff; }
  // Convert a 32 bit signed value from file's native byte order
  int Get32s(void * Long){
    return (motorolaOrder)?
      ((( char *)Long)[0] << 24) | (((uchar *)Long)[1] << 16) | (((uchar *)Long)[2] << 8 ) | (((uchar *)Long)[3] << 0 ) :
      ((( char *)Long)[3] << 24) | (((uchar *)Long)[2] << 16) | (((uchar *)Long)[1] << 8 ) | (((uchar *)Long)[0] << 0 ) ;
  }

  // Evaluate number, be it int, rational, or float from directory.
  double ConvertAnyFormat(void * ValuePtr, int Format){
    double Value = 0;

    switch(Format){
      case FMT_SBYTE:     Value = *(signed char *)ValuePtr;  break;
      case FMT_BYTE:      Value = *(uchar *)ValuePtr;        break;
      case FMT_USHORT:    Value = Get16u(ValuePtr);          break;
      case FMT_ULONG:     Value = Get32u(ValuePtr);          break;
      case FMT_URATIONAL:
      case FMT_SRATIONAL:{
        int Num = Get32s(ValuePtr);
        int Den = Get32s(4+(char *)ValuePtr);
        Value = (Den == 0)? 0 : double(Num/Den);
        break;
      }
      case FMT_SSHORT:    Value = (signed short)Get16u(ValuePtr);  break;
      case FMT_SLONG:     Value = Get32s(ValuePtr);                break;
      // Not sure if this is correct (never seen float used in Exif format)
      case FMT_SINGLE:    Value = (double)*(float *)ValuePtr;      break;
      case FMT_DOUBLE:    Value = *(double *)ValuePtr;             break;
      default:
        fprintf(stderr, "Illegal format code %d", Format);
    }

    return Value;
  }

  bool readJpegSections(FILE *fp){
    //parse the marker stream until SOS or EOI is seen
    if(fgetc(fp) != 0xff)  return false;
    if(fgetc(fp) != M_SOI) return false;

    uint status = 0;
    while(!status){
      Section sec;

      int marker = 0;

      for(int i=0 ; i<=16 ; i++){
        marker = fgetc(fp);
        if(marker) break;
      }
      if(!marker){
        fprintf(stderr,"too many padding bytes\n");
        status = 1;
        continue;
      }
/*
      for(int a=0;a<=16;a++){
        marker = fgetc(fp);
        if(marker != 0xff) break;
        if(a >= 16){
          fprintf(stderr, "too many padding bytes\n");
          status = 1;
          continue;
        }
      }
*/
      sec.type = marker;

      //read the length of the section
      int lh = fgetc(fp);
      int ll = fgetc(fp);
      int itemlen = (lh << 8) | ll;
      if(itemlen < 2){
        fprintf(stderr, "invalid marker\n");
        status = 1;
        continue;
      }
      sec.size = itemlen;

      sec.data = (uchar*)malloc(itemlen);
      //store first two pre-read bytes
      sec.data[0] = (uchar)lh;
      sec.data[1] = (uchar)ll;

      int got = fread(sec.data+2, 1, itemlen-2, fp);
      if(itemlen-2 != got){ //read the whole section
        fprintf(stderr, "Premature end of file?\n");
        status = 1;
        continue;
      }

      switch(marker){
        case M_SOS: //stop before hitting compressed data 
          status = 2;
          continue;
        case M_EOI: //in case it's a tables-only JPEG stream
          fprintf(stderr, "No image in jpeg!\n");
          status = 1;
          continue;
        case M_COM: //comment section
          //process_COM(data, itemlen);
          break;
        case M_JFIF:
          // Regular jpegs always have this tag, exif images have the exif
          // marker instead, althogh ACDsee will write images with both markers.
          // This program will re-create this marker on absence of exif marker,
          // hence no need to keep the copy from the file.
          free(sec.data);
          sec.data = NULL;
          break;
        case M_EXIF:
          // There can be different section using the same marker. Ignore all but "Exif" one
          if(memcmp(sec.data+2, "Exif", 4) == 0){
            if( !process_EXIF(sec.data, itemlen) ) status = 1;
          }else // Oterwise, discard this section
            free(sec.data);
            sec.data = NULL;
          break;
        case M_IPTC:
        case M_SOF0:
        case M_SOF1:
        case M_SOF2:
        case M_SOF3:
        case M_SOF5:
        case M_SOF6:
        case M_SOF7:
        case M_SOF9:
        case M_SOF10:
        case M_SOF11:
        case M_SOF13:
        case M_SOF14:
        case M_SOF15:
          process_SOFn(sec.data, marker);
          break;
        default: // Skip any other sections
          break;
      }
    }
    return (status==2);
  }

  /*
    Process a COM marker.
    we must guard against random junk and varying newline representations.
  */
/* UNUSED
  void process_COM(const uchar *data, uint length){
    char Comment[MAX_COMMENT_SIZE+1];
    int nch = 0;

    length = max(length, MAX_COMMENT_SIZE); //truncate if it won't fit in our structure

    for(int a=2 ; a<length ; a++){
      uchar ch = data[a];
      if( (ch=='\r') && (data[a+1]=='\n') ) continue; //remove cr followed by lf
      Comment[nch++] = ( (ch>=32) || (ch=='\n') || (ch=='\t') )? char(ch) : '?';
    }
    Comment[nch] = '\0'; //null terminate
  }
*/

  // Process exif format directory, as used by Cannon maker note
  void ProcessCanonMakerNoteDir(unsigned char * DirStart, unsigned char * OffsetBase, unsigned ExifLength){
    int NumDirEntries;

    NumDirEntries = Get16u(DirStart);
    #define DIR_ENTRY_ADDR(Start, Entry) (Start+2+12*(Entry))

    uchar *DirEnd;
    DirEnd = DIR_ENTRY_ADDR(DirStart, NumDirEntries);
    if(DirEnd > (OffsetBase+ExifLength)){
      fprintf(stderr, "Illegally sized exif makernote subdir (%d entries)", NumDirEntries);
      return;
    }

    for(int de=0 ; de<NumDirEntries ; de++){
      int Tag, Format, Components;
      uchar *ValuePtr;
      int ByteCount;
      uchar *DirEntry;
      DirEntry = DIR_ENTRY_ADDR(DirStart, de);

      Tag = Get16u(DirEntry);
      Format = Get16u(DirEntry+2);
      Components = Get32u(DirEntry+4);

      if((Format-1) >= NUM_FORMATS){
        // (-1) catches illegal zero case as unsigned underflows to positive large.
        fprintf(stderr, "Illegal number format %d for tag %04x", Format, Tag);
        continue;
      }
      if((unsigned)Components > 0x10000){
        fprintf(stderr, "Illegal number of components %d for tag %04x", Components, Tag);
        continue;
      }

      ByteCount = Components * BytesPerFormat[Format];

      if (ByteCount > 4){
        unsigned OffsetVal;
        OffsetVal = Get32u(DirEntry+8);
        // If its bigger than 4 bytes, the dir entry contains an offset.
        if (OffsetVal+ByteCount > ExifLength){
            // Bogus pointer offset and / or bytecount value
            fprintf(stderr, "Illegal value pointer for tag %04x", Tag);
            continue;
        }
        ValuePtr = OffsetBase+OffsetVal;
      }else{
          // 4 bytes or less and value is in the dir entry itself
          ValuePtr = DirEntry+8;
      }

      if( (Tag == 1) && (Components > 16) ){
        int IsoCode = Get16u(ValuePtr + 16*sizeof(unsigned short));
        if( (IsoCode >= 16) && (IsoCode <= 24) )
          ImageInfo.ISOequivalent = 50 << (IsoCode-16);
      }

      if( (Tag == 4) && (Format == FMT_USHORT) ){
        if(Components > 7){
          int WhiteBalance = Get16u(ValuePtr + 7*sizeof(unsigned short));
          switch(WhiteBalance){
            // 0=Auto, 6=Custom
            case 1: ImageInfo.LightSource = 1; break; // Sunny
            case 2: ImageInfo.LightSource = 1; break; // Cloudy
            case 3: ImageInfo.LightSource = 3; break; // Thungsten
            case 4: ImageInfo.LightSource = 2; break; // Fourescent
            case 5: ImageInfo.LightSource = 4; break; // Flash
          }
        }
        if( (Components > 19) && (ImageInfo.Distance <= 0) ){
          // Indicates the distance the autofocus camera is focused to.
          // Tends to be less accurate as distance increases.
          int temp_dist = Get16u(ValuePtr + 19*sizeof(unsigned short));
          ImageInfo.Distance = (temp_dist != 65535)? (float)temp_dist/100 : -1;
        }
      }
    }
  }

  // Process maker note - to the limited extent that its supported.
  void ProcessMakerNote(unsigned char * ValuePtr, unsigned char * OffsetBase, unsigned ExifLength){
    if(strstr(ImageInfo.CameraMake, "Canon"))
      ProcessCanonMakerNoteDir(ValuePtr, OffsetBase, ExifLength);
  }

  /*
    Process a SOFn marker. This is useful for the image dimensions
    JPEG image is data[7] color components, data[2] bits per sample
  */
  void process_SOFn(const uchar *data, int marker){
    //data[2] contains the data precision value
    ImageInfo.Height = Get16m(data+3);
    ImageInfo.Width = Get16m(data+5);
    int num_components = data[7];
    ImageInfo.IsColor = (num_components == 3);
    ImageInfo.Process = marker;
  }

  /*
    Process a EXIF marker.
    Describes all the drivel that most digital cameras include...
  */
  bool process_EXIF(uchar *ExifSection, uint length){
    FocalplaneXRes = 0;
    FocalplaneUnits = 0;
    ExifImageWidth = 0;
    NumOrientations = 0;

    // Check the EXIF header component
    static char *ExifHeader = "Exif\0\0";
    if(memcmp(ExifSection+2, ExifHeader,6)){
      fprintf(stderr, "Incorrect Exif header");
      return false;
    }
    if(memcmp(ExifSection+8,"II",2) == 0)      motorolaOrder = false; //Exif section in Intel order
    else if(memcmp(ExifSection+8,"MM",2) == 0) motorolaOrder = true;  //Exif section in Motorola order
    else{
      fprintf(stderr, "Invalid Exif alignment marker");
      return false;
    }
    // Check the next value for correctness
    if(Get16u(ExifSection+10) != 0x2a){
      fprintf(stderr, "Invalid Exif start (1)");
      return false;
    }

    const uint FirstOffset = Get32u(ExifSection+12);
    if(FirstOffset < 8 || FirstOffset > 16){
      // Usually set to 8, but other values valid too
      fprintf(stderr, "Suspicious offset of first IFD value");
      return false;
    }

    DirWithThumbnailPtrs = NULL;

    // First directory starts 16 bytes in. All offset are relative to 8 bytes in.
    if( !ProcessExifDir(ExifSection+8+FirstOffset, ExifSection+8, length-8, 0) ) return false;

    ImageInfo.ThumbnailAtEnd = (ImageInfo.ThumbnailOffset >= ImageInfo.LargestExifOffset);

    return true;
  }

  // Process one of the nested EXIF directories.
  bool ProcessExifDir(uchar *DirStart, uchar *OffsetBase, unsigned ExifLength, int NestingLevel){
    int NumDirEntries;
    unsigned ThumbnailOffset = 0;
    unsigned ThumbnailSize = 0;
    char IndentString[25];

    if(NestingLevel > 4){
      fprintf(stderr, "Maximum directory nesting exceeded (corrupt exif header)");
      return false;
    }

    memset(IndentString, ' ', 25);
    IndentString[NestingLevel * 4] = '\0';

    NumDirEntries = Get16u(DirStart);
#define DIR_ENTRY_ADDR(Start, Entry) (Start+2+12*(Entry))
    uchar *DirEnd = DIR_ENTRY_ADDR(DirStart, NumDirEntries);
    if(DirEnd+4 > (OffsetBase+ExifLength)){
      if(DirEnd+2 == OffsetBase+ExifLength || DirEnd == OffsetBase+ExifLength){
        // Version 1.3 of jhead would truncate a bit too much.
        // This also caught later on as well.
      }else{
        fprintf(stderr, "Illegally sized exif subdirectory (%d entries)", NumDirEntries);
        return false;
      }
    }

    for(int de=0 ; de<NumDirEntries ; de++){
      int Tag, Format, Components;
      uchar *DirEntry;
      DirEntry = DIR_ENTRY_ADDR(DirStart, de);

      Tag = Get16u(DirEntry);
      Format = Get16u(DirEntry+2);
      Components = Get32u(DirEntry+4);

      if((Format-1) >= NUM_FORMATS){
        // (-1) catches illegal zero case as unsigned underflows to positive large.
        fprintf(stderr, "Illegal number format %d for tag %04x", Format, Tag);
        continue;
      }

      if ((unsigned)Components > 0x10000){
        fprintf(stderr, "Illegal number of components %d for tag %04x", Components, Tag);
        continue;
      }

      int ByteCount = Components * BytesPerFormat[Format];

      uchar *ValuePtr;
      if(ByteCount > 4){
        unsigned OffsetVal;
        OffsetVal = Get32u(DirEntry+8);
        // If its bigger than 4 bytes, the dir entry contains an offset.
        if(OffsetVal+ByteCount > ExifLength){
          // Bogus pointer offset and / or bytecount value
          fprintf(stderr, "Illegal value pointer for tag %04x", Tag,0);
          continue;
        }
        ValuePtr = OffsetBase+OffsetVal;
        if(OffsetVal > ImageInfo.LargestExifOffset) ImageInfo.LargestExifOffset = OffsetVal;
      }else{
        // 4 bytes or less and value is in the dir entry itself
        ValuePtr = DirEntry+8;
      }

      if(Tag == TAG_MAKER_NOTE){
        ProcessMakerNote(ValuePtr, OffsetBase, ExifLength);
        continue;
      }

      // Extract useful components of tag
      switch(Tag){
        case TAG_MAKE:
          strncpy(ImageInfo.CameraMake, (char*)ValuePtr, ByteCount < 31 ? ByteCount : 31);
          break;
        case TAG_MODEL:
          strncpy(ImageInfo.CameraModel, (char*)ValuePtr, ByteCount < 39 ? ByteCount : 39);
          break;
        case TAG_DATETIME_ORIGINAL:
          // If we get a DATETIME_ORIGINAL, we use that one.
          strncpy(ImageInfo.DateTime, (char*)ValuePtr, 19);
          // Fallthru...
        case TAG_DATETIME_DIGITIZED:
        case TAG_DATETIME:
          // If we don't already have a DATETIME_ORIGINAL, use whatever time fields we may have.
          if(!isdigit(ImageInfo.DateTime[0])) strncpy(ImageInfo.DateTime, (char*)ValuePtr, 19);
          if(ImageInfo.numDateTimeTags >= MAX_DATE_COPIES)
            fprintf(stderr, "More than %d date fields!  This is nuts", MAX_DATE_COPIES, 0);
          else
            ImageInfo.DateTimeOffsets[ImageInfo.numDateTimeTags++] = (char*)ValuePtr - (char*)OffsetBase;
          break;
        case TAG_WINXP_COMMENT:
          // We already have a jpeg comment (probably windows comment), skip this one.
          if(ImageInfo.Comments[0]) break;
          if(ByteCount > 1){
            if(ByteCount > MAX_COMMENT_SIZE) ByteCount = MAX_COMMENT_SIZE;
            memcpy(ImageInfo.Comments, ValuePtr, ByteCount);
            ImageInfo.CommentWidchars = ByteCount/2;
          }
          break;
        case TAG_USERCOMMENT:{
          // We already have a jpeg comment (probably windows comment), skip this one.
          if(ImageInfo.Comments[0]) break;
          // Comment is often padded with trailing spaces. Remove these first.
          int a = ByteCount-1;
          while( a && ((ValuePtr)[a]==' ') ) (ValuePtr)[a--] = '\0';

          // Copy the comment
          if(memcmp(ValuePtr, "ASCII",5) == 0){
            for(a=5 ; a<10 ; a++){
              uchar c = (ValuePtr)[a];
              if( (c!='\0') && (c!=' ') ){
                strncpy(ImageInfo.Comments, (char *)ValuePtr+a, 199);
                break;
              }
            }
          }else{
            strncpy(ImageInfo.Comments, (char *)ValuePtr, MAX_COMMENT_SIZE-1);
          }
          break;
        }
        case TAG_FNUMBER:
          // Simplest way of expressing aperture, so I trust it the most.
          // (overwrite previously computd value if there is one)
          ImageInfo.ApertureFNumber = T(ConvertAnyFormat(ValuePtr, Format));
          break;
        case TAG_APERTURE:
        case TAG_MAXAPERTURE:
          // More relevant info always comes earlier, so only use this field if we don't 
          // have appropriate aperture information yet.
          if(ImageInfo.ApertureFNumber == 0)
            ImageInfo.ApertureFNumber = T(exp(ConvertAnyFormat(ValuePtr, Format)*log(2.0)*0.5));
          break;
        case TAG_FOCALLENGTH:
          // Nice digital cameras actually save the focal length as a function
          // of how farthey are zoomed in.
          ImageInfo.FocalLength = T(ConvertAnyFormat(ValuePtr, Format));
          break;
        case TAG_SUBJECT_DISTANCE:
          // Inidcates the distacne the autofocus camera is focused to.
          // Tends to be less accurate as distance increases.
          ImageInfo.Distance = T(ConvertAnyFormat(ValuePtr, Format));
          break;
        case TAG_EXPOSURETIME:
          // Simplest way of expressing exposure time, so I trust it most.
          // (overwrite previously computd value if there is one)
          ImageInfo.ExposureTime = T(ConvertAnyFormat(ValuePtr, Format));
          break;
        case TAG_SHUTTERSPEED:
          // More complicated way of expressing exposure time, so only use
          // this value if we don't already have it from somewhere else.
          if(ImageInfo.ExposureTime == 0)
            ImageInfo.ExposureTime = T(1/exp(ConvertAnyFormat(ValuePtr, Format)*log(2.0)));
          break;
        case TAG_FLASH:
          ImageInfo.FlashUsed=(int)ConvertAnyFormat(ValuePtr, Format);
          break;
        case TAG_ORIENTATION:
          if(NumOrientations >= 2){
            // Can have another orientation tag for the thumbnail, but if there's a third one, things are stringae.
            fprintf(stderr, "More than two orientation tags!");
            break;
          }
          OrientationPtr[NumOrientations] = ValuePtr;
          OrientationNumFormat[NumOrientations] = Format;
          if(NumOrientations == 0)
            ImageInfo.Orientation = (int)ConvertAnyFormat(ValuePtr, Format);
          if(ImageInfo.Orientation < 0 || ImageInfo.Orientation > 8){
            fprintf(stderr, "Undefined rotation value %d", ImageInfo.Orientation);
            ImageInfo.Orientation = 0;
          }
          NumOrientations += 1;
          break;
        case TAG_EXIF_IMAGELENGTH:
        case TAG_EXIF_IMAGEWIDTH:{
          // Use largest of height and width to deal with images that have been
          // rotated to portrait format.
          int a = (int)ConvertAnyFormat(ValuePtr, Format);
          if (ExifImageWidth < a) ExifImageWidth = a;
          break;
        }
        case TAG_FOCAL_PLANE_XRES:
          FocalplaneXRes = ConvertAnyFormat(ValuePtr, Format);
          break;
        case TAG_FOCAL_PLANE_UNITS:
          switch((int)ConvertAnyFormat(ValuePtr, Format)){
            case 1: FocalplaneUnits = 25.4; break; // inch
            case 2: 
              // According to the information I was using, 2 means meters.
              // But looking at the Cannon powershot's files, inches is the only
              // sensible value.
              FocalplaneUnits = 25.4;
              break;
            case 3: FocalplaneUnits = 10;   break;  // centimeter
            case 4: FocalplaneUnits = 1;    break;  // millimeter
            case 5: FocalplaneUnits = .001; break;  // micrometer
          }
          break;
        case TAG_EXPOSURE_BIAS:
          ImageInfo.ExposureBias = (float)ConvertAnyFormat(ValuePtr, Format);
          break;
        case TAG_WHITEBALANCE:
          ImageInfo.Whitebalance = (int)ConvertAnyFormat(ValuePtr, Format);
          break;
        case TAG_LIGHT_SOURCE:
          ImageInfo.LightSource = (int)ConvertAnyFormat(ValuePtr, Format);
          break;
        case TAG_METERING_MODE:
          ImageInfo.MeteringMode = (int)ConvertAnyFormat(ValuePtr, Format);
          break;
        case TAG_EXPOSURE_PROGRAM:
          ImageInfo.ExposureProgram = (int)ConvertAnyFormat(ValuePtr, Format);
          break;
        case TAG_EXPOSURE_INDEX:
          if(ImageInfo.ISOequivalent == 0){
            // Exposure index and ISO equivalent are often used interchangeably,
            // so we will do the same in jhead.
            // http://photography.about.com/library/glossary/bldef_ei.htm
            ImageInfo.ISOequivalent = (int)ConvertAnyFormat(ValuePtr, Format);
          }
          break;
        case TAG_EXPOSURE_MODE:
          ImageInfo.ExposureMode = (int)ConvertAnyFormat(ValuePtr, Format);
          break;
        case TAG_ISO_EQUIVALENT:
          ImageInfo.ISOequivalent = (int)ConvertAnyFormat(ValuePtr, Format);
          // Fixes strange encoding on some older digicams.
          if( ImageInfo.ISOequivalent < 50 ) ImageInfo.ISOequivalent *= 200;
          break;
        case TAG_DIGITALZOOMRATIO:
          ImageInfo.DigitalZoomRatio = (float)ConvertAnyFormat(ValuePtr, Format);
          break;
        case TAG_THUMBNAIL_OFFSET:
          ThumbnailOffset = (unsigned)ConvertAnyFormat(ValuePtr, Format);
          DirWithThumbnailPtrs = DirStart;
          break;
        case TAG_THUMBNAIL_LENGTH:
          ThumbnailSize = (unsigned)ConvertAnyFormat(ValuePtr, Format);
          ImageInfo.ThumbnailSizeOffset = ValuePtr-OffsetBase;
          break;
        case TAG_EXIF_OFFSET:
        case TAG_INTEROP_OFFSET:{
          uchar *SubdirStart = OffsetBase + Get32u(ValuePtr);
          if(SubdirStart < OffsetBase || SubdirStart > OffsetBase+ExifLength)
            fprintf(stderr, "Illegal exif or interop offset directory link");
          else
            ProcessExifDir(SubdirStart, OffsetBase, ExifLength, NestingLevel+1);
          break;
        }
        case TAG_GPSINFO:{
          uchar *SubdirStart = OffsetBase + Get32u(ValuePtr); 
          if(SubdirStart < OffsetBase || SubdirStart > OffsetBase+ExifLength)
            fprintf(stderr, "Illegal GPS directory link");
          else
            ;//ProcessGpsInfo(SubdirStart, ByteCount, OffsetBase, ExifLength);
          break;
        }
        case TAG_FOCALLENGTH_35MM:
          // The focal length equivalent 35 mm is a 2.2 tag (defined as of April 2002)
          // if its present, use it to compute equivalent focal length instead of 
          // computing it from sensor geometry and actual focal length.
          ImageInfo.FocalLength35mmEquiv = (unsigned)ConvertAnyFormat(ValuePtr, Format);
          break;
        case TAG_DISTANCE_RANGE:
          // Three possible standard values:
          //   1 = macro, 2 = close, 3 = distant
          ImageInfo.DistanceRange = (int)ConvertAnyFormat(ValuePtr, Format);
          break;
      }
    }

    // In addition to linking to subdirectories via exif tags, 
    // there's also a potential link to another directory at the end of each
    // directory.  this has got to be the result of a committee!
    unsigned char * SubdirStart;
    unsigned Offset;

    if(DIR_ENTRY_ADDR(DirStart, NumDirEntries) + 4 <= OffsetBase+ExifLength){
      Offset = Get32u(DirStart+2+12*NumDirEntries);
      if (Offset){
        SubdirStart = OffsetBase + Offset;
        if(SubdirStart > OffsetBase+ExifLength || SubdirStart < OffsetBase){
          if (SubdirStart > OffsetBase && SubdirStart < OffsetBase+ExifLength+20){
            // Jhead 1.3 or earlier would crop the whole directory!
            // As Jhead produces this form of format incorrectness, 
            // I'll just let it pass silently
          }else{
            fprintf(stderr, "Illegal subdirectory link");
          }
        }else{
          if(SubdirStart <= OffsetBase+ExifLength)
            ProcessExifDir(SubdirStart, OffsetBase, ExifLength, NestingLevel+1);
        }
        if(Offset > ImageInfo.LargestExifOffset){
          ImageInfo.LargestExifOffset = Offset;
        }
      }
    }else{
      // The exif header ends before the last next directory pointer.
    }

    if(ThumbnailOffset){
      ImageInfo.ThumbnailAtEnd = false;
      if(ThumbnailOffset <= ExifLength){
        if(ThumbnailSize > ExifLength-ThumbnailOffset){
          // If thumbnail extends past exif header, only save the part that
          // actually exists.  Canon's EOS viewer utility will do this - the
          // thumbnail extracts ok with this hack.
          ThumbnailSize = ExifLength-ThumbnailOffset;
        }
        // The thumbnail pointer appears to be valid.  Store it.
        ImageInfo.ThumbnailOffset = ThumbnailOffset;
        ImageInfo.ThumbnailSize = ThumbnailSize;
      }
    }
  }

};

#endif