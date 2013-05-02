#ifndef GETFILEVALUE
#define GETFILEVALUE

bool getFileValue( const char* inFile, const char* parameterName, float& returnValue );

bool getFileValueWithError( const char* inFile, const char* parameterName, float& returnValue, float& returnError );

#endif
