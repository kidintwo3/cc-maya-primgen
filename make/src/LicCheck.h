#ifndef LICCHECK_H
#define LICCHECK_H

#include <maya/MString.h>
#include <maya/MGlobal.h>

#include <string>
#include <vector>

#include "curl/curl.h"

std::string buffer;


size_t curl_write( void *ptr, size_t size, size_t nmemb, void *stream)
{
    buffer.append((char*)ptr, size*nmemb);
    return size*nmemb;
}


bool checkLic(MString license_key)
{

	int result = false;
	
	CURL *curl = curl_easy_init();
	CURLcode res;
	std::string readBuffer;

	MString product_permalink = "primGen";
	// MString license_key = "D0F9F4D8-50FD4FE3-A411C173-FB6C77C7";

	if(curl) 
	{
		curl_easy_setopt(curl, CURLOPT_URL, "https://api.gumroad.com/v2/licenses/verify");

		curl_easy_setopt(curl, CURLOPT_POSTFIELDS, MString("product_permalink=" + product_permalink + "&license_key=" + license_key).asChar() );

		curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curl_write);

		/* Perform the request, res will get the return code */ 
		res = curl_easy_perform(curl);


		const char* primitiveStrings = "success\":true";


		int str_result = MString(buffer.c_str()).indexW(primitiveStrings);

		if (str_result != -1)
		{
			result = true;
		}
		

		/* always cleanup */ 
		curl_easy_cleanup(curl);

	}

	curl_global_cleanup();


	return result;

}


#endif
