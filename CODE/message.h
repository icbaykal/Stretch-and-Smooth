#ifndef ICB_MESSAGE_H
#define ICB_MESSAGE_H


#include <windows.h>
#include <string>

using namespace std;

//------------MESSAGE ----------------
#define MESG_BUF_SIZE	2048
class Message
{
	char messagebuf[MESG_BUF_SIZE];
	wchar_t wchbuf[MESG_BUF_SIZE];
	int result;
	int fprecision;
	int FontCol;
	HWND Editwin;
public:
	virtual ~Message(){}
	Message();
	void SetWindow(HWND c);
	void SetWindowColor(int BGR);
	void SetFontColor(int BGR);
	int Post(const char *ch);
	int Post(wchar_t *wch);
	int Post(char *ch1,char *ch2);
	int Post(int i);
	int Post(const char *ch,int i);
	int Post(const char *ch,int i,int j);
	int Post(const char *ch,int i,int j,int k);
	int Post(unsigned int i);
	int Post(const char *ch,unsigned int i);
	int Post(const char *ch,unsigned int i,unsigned int j);
	int Post(const char *ch,int unsigned i,int unsigned j,int unsigned k);
	int Post(float f);
	int Post(const char *ch,float f);
	int Post(const char *ch,float f,float g);
	int Post(double d);
	int Post(const char *ch,double d);
	int Post(const char *ch,double d,double e);
	int Post(const char *ch,char *d);
	int Post(string &s);
	int Post(const char *ch,int i,float f);
	int Post(const char *ch,float f,int i);
	int Post(const char *ch,int i,double d);
	int Post(const char *ch,double d,int i);
	//int Post(...); will work on this someday
	void Clear();
};

#endif
