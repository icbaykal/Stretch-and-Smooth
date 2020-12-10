#include <windows.h>
#include <string>
#include "message.h"
Message msg;

void Message::SetWindow(HWND c){Editwin=c;return;}
Message::Message(){messagebuf[0]=0;fprecision=4;FontCol=0;}
void Message::Clear(){
	memset(messagebuf,0,2048);SendMessage (Editwin, EM_SETSEL, (WPARAM)0, 
		(LPARAM)-1);SendMessage (Editwin, EM_REPLACESEL, 0, (LPARAM)  messagebuf);	
	//SetFontColor(FontCol);
UpdateWindow(Editwin);}//SetWindowText(Editwin,NULL);}

int Message::Post(const char *text)
{
    int len=strlen(text)+1;
	unsigned short *wc=new unsigned short[len+1];
	if(wc==NULL)return -1;
	//the length odf the text inside the window
	int textlen=GetWindowTextLength(Editwin);
	MultiByteToWideChar(CP_ACP, 0, (LPCSTR)text, len, (LPWSTR)wc, len);
	wc[len]=0;
	SendMessageW(Editwin, EM_SETSEL, (WPARAM)textlen, (LPARAM)textlen);
	SendMessageW(Editwin, EM_REPLACESEL, 0, (LPARAM)  wc);
	SendMessageW(Editwin,WM_VSCROLL,SB_PAGEDOWN,0);
	delete []wc;
	return(0);
}
int Message::Post(wchar_t *wch)
{    
    //the length of the text inside the window
	int textlen=GetWindowTextLength(Editwin);
	//select end of the text in the window
	SendMessageW(Editwin, EM_SETSEL, (WPARAM)textlen, (LPARAM)textlen);
	//replae non existing selection with the new text or append it.
	SendMessageW(Editwin, EM_REPLACESEL, 0, (LPARAM)  wch);
    return(0);    
}
//-----------CHAR -----------------
int Message::Post(char *ch1,char *ch2)
{
	messagebuf[0]=0; 
   sprintf(messagebuf,ch1);strcat(messagebuf,ch2); Post(messagebuf);
   return(1);
}
// ------------- INT  ----------------
int Message::Post(int i)
{
    wchbuf[0]=0;swprintf(wchbuf,L"%d",i);Post(wchbuf);
    return(0);
}
int Message::Post(const char *ch,int i)
{
   messagebuf[0]=0; 
   sprintf(messagebuf,ch,i); Post(messagebuf);
   return(1);
}
int Message::Post(const char *ch,int i,int j)
{	
	messagebuf[0]=0; sprintf(messagebuf,ch,i,j); Post(messagebuf);
	return(1);
}
int Message::Post(const char *ch,int i,int j,int k)
{	
	messagebuf[0]=0; sprintf(messagebuf,ch,i,j,k); Post(messagebuf);
	return(1);
}
//------------------UNSIGNED INT -------------------
int Message::Post(unsigned int i)
{
	sprintf(messagebuf,"%u",i); Post(messagebuf);
	return(1);
}
int Message::Post(const char *ch,unsigned int i)
{
	messagebuf[0]=0; sprintf(messagebuf,ch,i); Post(messagebuf);
	return(1);
}
int Message::Post(const char *ch,unsigned int i,unsigned int j)
{	
	messagebuf[0]=0; sprintf(messagebuf,ch,i,j); Post(messagebuf);
	return(1);
}
int Message::Post(const char *ch,unsigned int i,unsigned int j,unsigned int k)
{	
	messagebuf[0]=0; sprintf(messagebuf,ch,i,j,k); Post(messagebuf);
	return(1);
}
//----------------- FLOAT -----------------
int Message::Post(float f)
{
	sprintf(messagebuf,"%f",f); Post(messagebuf);
	return(1);
}
int Message::Post(const char *ch,float f)
{	
	messagebuf[0]=0; sprintf(messagebuf,ch,f); Post(messagebuf);
	return(1);
}
int Message::Post(const char *ch,float f,float g)
{	
	messagebuf[0]=0; sprintf(messagebuf,ch,f,g); Post(messagebuf);
	return(1);
}
//--------------------DOUBLE -----------------
int Message::Post(double d)
{
	sprintf(messagebuf,"%e",d); Post(messagebuf);
	return(1);
}
int Message::Post(const char *ch,double d)
{	
	messagebuf[0]=0; sprintf(messagebuf,ch,d); Post(messagebuf);
	return(1);
}
int Message::Post(const char *ch,double d,double e)
{	
	messagebuf[0]=0; sprintf(messagebuf,ch,d,e); Post(messagebuf);
	return(1);
}
//-----------------String --------------------
int Message::Post(const char *ch,char *d)
{	
	messagebuf[0]=0; sprintf(messagebuf,ch,d); Post(messagebuf);
	return(1);
}

int Message::Post(string &s)
{	
	const char *k=s.data();
	Post(k);
	return(1);
}
//-----------------INT % FLOAT --------------------
int Message::Post(const char *ch,int i,float f)
{	
	messagebuf[0]=0; sprintf(messagebuf,ch,i,f); Post(messagebuf);
	return(1);
}
int Message::Post(const char *ch,float f,int i)
{	
	messagebuf[0]=0; sprintf(messagebuf,ch,f,i); Post(messagebuf);
	return(1);
}
//-----------------INT % DOUBLE --------------------
int Message::Post(const char *ch,int i,double d)
{	
	messagebuf[0]=0; sprintf(messagebuf,ch,i,d); Post(messagebuf);
	return(1);
}
int Message::Post(const char *ch,double d,int i)
{	
	messagebuf[0]=0; sprintf(messagebuf,ch,d,i); Post(messagebuf);
	return(1);
}
