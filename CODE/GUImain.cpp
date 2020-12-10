#include <windows.h>  
#include <commctrl.h>
#include <chrono>
#include "message.h"
#include "kht.hpp"

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp> 
#include <opencv2/opencv.hpp>

using namespace cv;
extern Message msg;  //the message object

#define IDI_ICON1 101
#define IDC_OPEN				3000
#define IDC_BILINEAR			3001
#define IDC_CANNY				3500
#define IDC_CANNY_x3			3501
#define IDC_THIN				4000
#define IDC_STRETCH				4001
#define IDC_SMOOTH				4002
#define IDC_STD_HOUGH			4998
#define IDC_STD_KHT				4999
#define IDC_PRO_HOUGH			5000
#define IDC_CRCL_HOUGH			5001
#define IDC_STD_HOUGH_BIL		5002
#define IDC_KHT_BIL				5003
#define IDC_PRO_HOUGH_BIL		5004
#define IDC_CRCL_HOUGH_BIL		5005
#define IDC_STD_HOUGH_STRTCH	5013
#define IDC_KHT_STRTCH			5014
#define IDC_PRO_HOUGH_STRTCH	5015
#define IDC_CRCL_HOUGH_STRTCH	5016

#define IDC_BILIN_STD_HOUGH_TEST			5998
#define IDC_STRETCH_STD_HOUGH_TEST			5999
#define IDC_BILIN_KHT_TEST					6000
#define IDC_STRETCH_KHT_TEST				6001
#define IDC_BILIN_PROB_HOUGH_TEST			6002
#define IDC_STRETCH_PROB_HOUGH_TEST			6003
#define IDC_BILIN_CIRC_HOUGH_TEST			6004
#define IDC_STRETCH_CIRC_HOUGH_TEST			6005

#define IDC_HELP				7000
#define IDC_ABOUT				7001

std::double_t const DEGREES_TO_RADIANS = std::atan(1.0) / 45.0;
HINSTANCE hInst;   // current instance
HWND hWnd;    //parent window
HWND mainbmp, linesbmp, edgebmp, mesgbox;
HWND realparent, edgeparent;
HWND ProbHoughstatic, ProbHoughEd1, ProbHoughEd2, ProbHoughEd3, ProbHoughEd4;
HWND StdHoughstatic, StdHoughEd1, StdHoughEd2, StdHoughEd3, StdHoughEd4;
HDC hDC;
HGLRC hRC;

LPCTSTR lpszAppName = "Probabilistic Hough";
LPCTSTR lpszTitle = "PROBABILISTIC HOUGH";
LRESULT CALLBACK WndProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam);

char mes[64], tmes[32];
Mat image, image_x3,imcolor;
Mat copyofimg;
Mat edges, specialthin, edges_x3;
Mat stretch;
Mat smooth;
int blinearflag = 0;

class STMImage
{
public:
	HBITMAP HBMRGB;
	BITMAPINFO info;
	unsigned char *buf;
	unsigned char *RGB;
	unsigned char *pieces;
	int *seg;
	int width;
	int height;
	unsigned char *piecetable;
	int noofsegments;
	int validsegments;
	STMImage() { buf = NULL; RGB = NULL; pieces = NULL; piecetable = NULL; }
	void freebuffers() { free(buf); free(RGB); free(pieces); free(piecetable); free(seg); }
	~STMImage() { freebuffers(); }
};
STMImage i;

int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
	LPTSTR lpCmdLine, int nCmdShow)
{
	MSG      mesg;
	WNDCLASSEX wc;
	HMENU MainMenu, FileMenu,OrigMetMenu,KHTMenu, MorphMenu,TransMenu, SpeedMenu, HelpMenu;
	MainMenu = CreateMenu();	
	FileMenu = CreatePopupMenu();
	OrigMetMenu=CreatePopupMenu();
	MorphMenu = CreatePopupMenu();
	TransMenu = CreatePopupMenu();
	SpeedMenu = CreatePopupMenu();
	HelpMenu = CreatePopupMenu();
	//FILE Menu ------------
	AppendMenu(FileMenu, MF_STRING, IDC_OPEN, "Open");		
	AppendMenu(MainMenu, MF_POPUP, (UINT_PTR)FileMenu, "FILE");
	//ORIGINAL METHOD Menu----------------
	AppendMenu(OrigMetMenu, MF_STRING, IDC_CANNY, "Canny Edge");
	AppendMenu(OrigMetMenu, MF_STRING, IDC_STD_HOUGH, "Standard Hough Tr.");
	AppendMenu(OrigMetMenu, MF_STRING, IDC_STD_KHT, "Kernel Based H.Tr.");
	AppendMenu(OrigMetMenu, MF_STRING, IDC_PRO_HOUGH, "Prob. Hough Tr.");
	AppendMenu(OrigMetMenu, MF_STRING, IDC_CRCL_HOUGH, "Circular Hough Tr.");
	AppendMenu(MainMenu, MF_POPUP, (UINT_PTR)OrigMetMenu, "ORIGINAL METHOD");
	// BILINEAR TRANSFORM Menu ------------
	AppendMenu(TransMenu, MF_STRING, IDC_BILINEAR, "Bilinear x3");
	AppendMenu(TransMenu, MF_STRING, IDC_CANNY_x3, "Canny Edge");
	AppendMenu(TransMenu, MF_STRING, IDC_STD_HOUGH_BIL, "Standard Hough Tr.");
	AppendMenu(TransMenu, MF_STRING, IDC_KHT_BIL, "Kernel Based H.Tr.");
	AppendMenu(TransMenu, MF_STRING, IDC_PRO_HOUGH_BIL, "Prob. Hough Tr.");
	AppendMenu(TransMenu, MF_STRING, IDC_CRCL_HOUGH_BIL, "Circular Hough Tr.");
	AppendMenu(MainMenu, MF_POPUP, (UINT_PTR)TransMenu, "BILINEAR TRANSFORM");
	// MORPHOLOGIC STRETCH Menu ------------------
	AppendMenu(MorphMenu, MF_STRING, IDC_CANNY, "Canny Edge");
	AppendMenu(MorphMenu, MF_STRING, IDC_THIN, "Special Thin");
	AppendMenu(MorphMenu, MF_STRING, IDC_STRETCH, "Stretch");
	AppendMenu(MorphMenu, MF_STRING, IDC_SMOOTH, "Smooth");
	AppendMenu(MorphMenu, MF_STRING, IDC_STD_HOUGH_STRTCH, "Standard Hough Tr.");
	AppendMenu(MorphMenu, MF_STRING, IDC_KHT_STRTCH, "Kernel Based H.Tr.");
	AppendMenu(MorphMenu, MF_STRING, IDC_PRO_HOUGH_STRTCH, "Prob. Hough Tr.");
	AppendMenu(MorphMenu, MF_STRING, IDC_CRCL_HOUGH_STRTCH, "Circular Hough Tr.");
	AppendMenu(MainMenu, MF_POPUP, (UINT_PTR)MorphMenu, "MORPHOLOGIC STRETCH");
	//SPEED TESTS Menu -----------------------	
	AppendMenu(SpeedMenu, MF_STRING, IDC_BILIN_STD_HOUGH_TEST, "Bilinear Stretch Stand. Hough Speed Test");
	AppendMenu(SpeedMenu, MF_STRING, IDC_STRETCH_STD_HOUGH_TEST, "Morphologic Stretch Stand. Hough Speed Test");
	AppendMenu(SpeedMenu, MF_STRING, IDC_BILIN_KHT_TEST, "Bilinear Stretch KHT Speed Test");
	AppendMenu(SpeedMenu, MF_STRING, IDC_STRETCH_KHT_TEST, "Morphologic Stretch KHT Speed Test");
	AppendMenu(SpeedMenu, MF_STRING, IDC_BILIN_PROB_HOUGH_TEST, "Bilinear Stretch Prob. Hough Speed Test");
	AppendMenu(SpeedMenu, MF_STRING, IDC_STRETCH_PROB_HOUGH_TEST, "Morphologic Stretch Prob. Hough Speed Test");
	AppendMenu(SpeedMenu, MF_STRING, IDC_BILIN_CIRC_HOUGH_TEST, "Bilinear Stretch Circular Hough Speed Test");
	AppendMenu(SpeedMenu, MF_STRING, IDC_STRETCH_CIRC_HOUGH_TEST, "Morphologic Stretch Circular Hough Speed Test");
	//AppendMenu(SpeedMenu, MF_STRING, IDC_SMOOTHTEST, "Smooth");
	AppendMenu(MainMenu, MF_POPUP, (UINT_PTR)SpeedMenu, "SPEED TESTS");
	AppendMenu(HelpMenu, MF_STRING, IDC_HELP, "Help");
	AppendMenu(HelpMenu, MF_STRING, IDC_ABOUT, "About");
	AppendMenu(MainMenu, MF_POPUP, (UINT_PTR)HelpMenu,"HELP");



	wc.style = CS_HREDRAW | CS_VREDRAW;
	wc.lpfnWndProc = (WNDPROC)WndProc;
	wc.cbClsExtra = 0;
	wc.cbWndExtra = 0;
	wc.hInstance = hInstance;
	wc.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_ICON1));
	wc.hCursor = LoadCursor(NULL, IDC_ARROW);
	wc.hbrBackground = (HBRUSH)(COLOR_WINDOW);
	wc.lpszMenuName = lpszAppName;
	wc.lpszClassName = lpszAppName;
	wc.cbSize = sizeof(WNDCLASSEX);
	wc.hIconSm = (HICON)LoadImage(hInstance, lpszAppName,
		IMAGE_ICON, 16, 16,
		LR_DEFAULTCOLOR);

	if (!RegisterClassEx(&wc))
		return(FALSE);

	hInst = hInstance;
	hWnd = CreateWindowEx(0, lpszAppName,
		lpszTitle,
		WS_OVERLAPPEDWINDOW,
		60, 0,
		1885, 1050,
		NULL,
		MainMenu,
		hInstance,
		NULL
	);
	//Initialize Common Controls  




	if (!hWnd)
		return(FALSE);
	//---------------------------------
	mainbmp = CreateWindow(WC_STATIC, "", WS_CHILD | WS_VISIBLE | SS_BITMAP,
		3, 0, 1200, 800, hWnd, NULL, hInst, NULL);
	cv::namedWindow("My Window", 1);//,CV_WINDOW_AUTOSIZE);
	HWND hWnd2 = (HWND)cvGetWindowHandle("My Window");
	realparent = GetParent(hWnd2);
	::SetParent(hWnd2, mainbmp);
	Mat temp;
	temp = Mat::zeros(480, 480, CV_8U);
	cv::imshow("My Window", temp);
	ShowWindow(realparent, SW_HIDE);
	temp = Mat::zeros(1, 1, CV_8U);
	/*edgebmp = CreateWindow(WC_STATIC, "", WS_CHILD | WS_VISIBLE | SS_BITMAP,
		10, 515, 640, 520, hWnd, NULL, hInst, NULL);
	cv::namedWindow("Edge Window", 1);//,CV_WINDOW_AUTOSIZE);
	HWND hWnd3 = (HWND)cvGetWindowHandle("Edge Window");
	edgeparent = GetParent(hWnd3);
	::SetParent(hWnd3, edgebmp);*/
	//-----------------------------

	ShowWindow(hWnd, nCmdShow);
	UpdateWindow(hWnd);

	while (GetMessage(&mesg, NULL, 0, 0))
	{
		TranslateMessage(&mesg);
		DispatchMessage(&mesg);
	}

	return(mesg.wParam);
}


int FillBMPInf(BITMAPINFO *info, int x, int y, int bits)
{
	info->bmiHeader.biBitCount = bits;
	info->bmiHeader.biClrImportant = 0;
	info->bmiHeader.biClrUsed = 0;
	info->bmiHeader.biCompression = BI_RGB;
	info->bmiHeader.biHeight = -y;
	info->bmiHeader.biPlanes = 1;
	info->bmiHeader.biSize = sizeof(BITMAPINFO);
	info->bmiHeader.biSizeImage = (x*y*bits) / 8;
	if (info->bmiHeader.biSizeImage < 0)
		info->bmiHeader.biSizeImage = -info->bmiHeader.biSizeImage;
	info->bmiHeader.biWidth = x;
	info->bmiHeader.biXPelsPerMeter = 0;
	info->bmiHeader.biYPelsPerMeter = 0;
	return(0);
}


int Line(STMImage &i, int xstart, int ystart, int xend, int yend, int color)
{
	unsigned int *RGB = (unsigned int *)i.RGB;
	double slope, dy;
	int x, y, temp;
	
	if (abs(xend - xstart) > abs(yend - ystart))
	{

		slope = (double)(yend - ystart);
		slope = slope / (double)(xend - xstart);
		dy = 0.0;
		y = ystart;
		for (x = xstart;x <= xend;x++)
		{
			if (y<2 || y>i.height-4) continue;
			if (x<2 || x>i.width - 4) continue;
			RGB[y*i.width + x] = color;
			RGB[(y + 1)*i.width + x] = color;
			RGB[(y + 2)*i.width + x] = color;
			dy += ((double)x*slope + 0.5);
			y = ystart + ((double)(x - xstart)*slope + 0.5);
		}
	}
	else
	{

		slope = (double)(xend - xstart) / (double)(yend - ystart);
		if (yend < ystart)
		{
			temp = yend;
			yend = ystart;
			ystart = temp;
			temp = xend;
			xend = xstart;
			xstart = temp;
		}
		for (y = ystart;y <= yend;y++)
		{
			x = xstart + ((double)(y - ystart)*slope + 0.5);
			if (y<2 || y>i.height - 4) continue;
			if (x<2 || x>i.width-4) continue;
			RGB[y*i.width + x] = color;
			RGB[y*i.width + x + 1] = color;
			RGB[y*i.width + x + 2] = color;
		}
	}
	return(0);
}

int Line(STMImage& i, double rho, double theta, int color)
{
	unsigned int* RGB = (unsigned int*)i.RGB;
	double slope, dy;
	int x, y, temp;
	//yatay
	if (abs(sin(theta))>0.5)
	{		
		
		for (x = -i.width; x <= i.width; x++)
		{
			//y = ystart +  ((double)x*slope + 0.5);
			y = (rho - (double)x * cos(theta)) / sin(theta);
			if (x > 0 && x < i.width && y>0 && y < (i.height - 1))
			{
				RGB[y * i.width + x] = color;
				RGB[(y + 1) * i.width + x] = color;
			}
			
		}
	}
	else
	{		
		for (y = -i.height; y <= i.height; y++)
		{
			x = (rho - (double)y * sin(theta)) / cos(theta);
			if (x > 0 && x < (i.width-1) && y>0 && y < i.height)
			{
				RGB[y * i.width + x] = color;
				RGB[y * i.width + x+1] = color;
			}
		}
	}
	return(0);
}

int Circle(STMImage& i, int xstart, int ystart, int r)
{
	double x;
	int xx, y;
	unsigned int* RGB = (unsigned int*)i.RGB;
	for (x = -r; x < -r + 3; x += 0.03)
	{
		y = sqrt(double(r * r - x * x)) + ystart;
		xx = x + xstart;
		if (y > 0 && y < i.height)
			if (xx > 0 && xx < i.width)
			{
				RGB[((int)y) * i.width + (int)x + xstart] = 0xffffff;
				RGB[y * i.width + (int)x + xstart + 1] = 0xffffff;

			}
		y = ystart * 2 - y;
		if (y < 1 || y >= i.height)continue;
		RGB[y * i.width + (int)x + xstart] = 0xffffff;
		RGB[y * i.width + (int)x + xstart + 1] = 0xffffff;

	}
	for (x = r - 3; x < r + 0.3; x += 0.03)
	{
		y = sqrt(double(r * r - x * x)) + ystart;
		xx = x + xstart;
		if (y > 0 && y < i.height)
			if (xx > 0 && xx < i.width)
			{
				RGB[((int)y) * i.width + (int)x + xstart] = 0xffffff;
				RGB[y * i.width + (int)x + xstart - 1] = 0xffffff;

			}
		y = ystart * 2 - y;
		if (y < 1 || y >= i.height)continue;
		RGB[y * i.width + (int)x + xstart] = 0xffffff;
		RGB[y * i.width + (int)x + xstart - 1] = 0xffffff;
	}
	//son nokta
	if (ystart > 0 && ystart < i.height && (xstart + r) >0 && (xstart + r) < i.width)
	{
		RGB[ystart * i.width + xstart + r] = 0xffffff;
		RGB[ystart * i.width + xstart + r - 1] = 0xffffff;
	}
	for (x = -r + 3 - 0.1; x < r - 3 + 0.1; x += 0.03)
	{
		y = sqrt(double(r * r - x * x)) + ystart;
		xx = x + xstart;
		if (y > 0 && y < i.height)
			if (xx > 0 && xx < i.width)
			{
				RGB[((int)y) * i.width + (int)x + xstart] = 0xffffff;
				RGB[y * i.width - i.width + (int)x + xstart] = 0xffffff;
			}
		y = ystart * 2 - y;
		if (y < 1 || y >= i.height)continue;
		RGB[y * i.width + (int)x + xstart] = 0xffffff;
		RGB[y * i.width + i.width + (int)x + xstart] = 0xffffff;
	}
	return(0);
}

void CopyToSTM(Mat &m, STMImage &i)
{
	for (int y = 0; y < i.height; y++)
		for (int x = 0; x < i.width; x++)
		{
			i.RGB[y * i.width * 4 + x * 4 + 1] = m.at<unsigned char>(y, x)*0.45;
			i.RGB[y*i.width * 4 + x * 4+2] = m.at<unsigned char>(y, x);
			i.RGB[y*i.width * 4 + x * 4 ] = m.at<unsigned char>(y, x);
			//u= edges.at<unsigned char>(x, y, 1);
		}
}

void ColorCopy(Mat& m, STMImage& i)
{
	for (int y = 0; y < i.height; y++)
		for (int x = 0; x < i.width; x++)
		{
			i.RGB[y * i.width * 4 + x * 4 + 1] = m.at<unsigned char>(y, x*3+1);
			i.RGB[y * i.width * 4 + x * 4 + 2] = m.at<unsigned char>(y, x*3+2);
			i.RGB[y * i.width * 4 + x * 4] = m.at<unsigned char>(y, x*3);
			//u= edges.at<unsigned char>(x, y, 1);
		}
}

void GetHoughParam(double &theta, int &thresh, double &minLinLen, double &maxLinLen)
{
	char ct[10];
	GetWindowTextA(ProbHoughEd1, ct, 10);
	 theta = atof(ct);
	GetWindowTextA(ProbHoughEd2, ct, 10);
	 thresh = atoi(ct);
	GetWindowTextA(ProbHoughEd3, ct, 10);
	 minLinLen = atof(ct);
	GetWindowTextA(ProbHoughEd4, ct, 10);
	 maxLinLen = atof(ct);
}
int MyStdHough(Mat& edges, vector <Vec2f>& lines)
{
	double theta, minLinLen, maxLinLen;
	int thresh;
	GetHoughParam(theta, thresh, minLinLen, maxLinLen);
	CopyToSTM(copyofimg, i);
	if (edges.rows < 3) return 0;
	HoughLines(edges, lines, 1.0, CV_PI / 180, 150, 0, 0);
	//HoughLines(edges, lines, 1.0, CV_PI / theta, thresh, minLinLen, maxLinLen);
	return lines.size();
}
int MyStdHoughx3(Mat& edges, vector <Vec2f>& lines)
{
	double theta, minLinLen, maxLinLen;
	int thresh;
	GetHoughParam(theta, thresh, minLinLen, maxLinLen);
	CopyToSTM(copyofimg, i);
	if (edges.rows < 3) return 0;
	HoughLines(edges, lines, 1.0, CV_PI / 180, 200, 0, 0);
	//HoughLines(edges, lines, 1.0, CV_PI / theta, thresh, minLinLen, maxLinLen);
	return lines.size();
}
int MyHough(Mat &edges, vector<Vec4i> &lines)
{
	double theta, minLinLen, maxLinLen;
	int thresh;
	GetHoughParam(theta, thresh, minLinLen, maxLinLen);
	CopyToSTM(copyofimg, i);
	
	if (edges.rows < 3) return 0;
	HoughLinesP(edges, lines, 1.0, CV_PI / theta, thresh, minLinLen, maxLinLen);	
	return lines.size();
}
int MyHoughx3(Mat &edges, vector<Vec4i> &lines)
{
	double theta, minLinLen, maxLinLen;
	int thresh;
	GetHoughParam(theta, thresh, minLinLen, maxLinLen);
	CopyToSTM(copyofimg, i);

	if (edges.rows < 3) return 0;
	//HoughLinesP(edges, lines, 1.0, CV_PI / theta, thresh, minLinLen, maxLinLen);
	HoughLinesP(edges, lines, 1.0, CV_PI / theta, thresh*1.25, minLinLen*1.25, maxLinLen*3.0);
	return lines.size();
}
void DrawLines(STMImage &img, vector<Vec4i> &lines)
{
	for (size_t j = 0; j < lines.size(); j++)
	{
		Vec4i l = lines[j];
		Line(img, l[0], l[1], l[2], l[3], 0xffffffff);
	}
 }

void DrawLines(STMImage& img, vector <Vec2f> lines)
{
	double rho,theta;
	Point pt1, pt2;
	double a, b;
	for (size_t ii = 0; ii < lines.size(); ii++)
	{
		rho = lines[ii][0]; theta = lines[ii][1];		
		
		Line(img, rho,theta, 0xffffffff);
	}
}

void DrawLines(STMImage& img, kht::ListOfLines lines, int maxlines)
{
	cv::Point p1, p2;
	int width = edges.cols, height = edges.rows;

	for (std::size_t j = 0; j != maxlines; ++j) {
		auto const& line = lines[j];
		std::double_t rho = line.rho;
		std::double_t theta = line.theta * DEGREES_TO_RADIANS;
		std::double_t cos_theta = cos(theta), sin_theta = sin(theta);

		// Convert from KHT to OpenCV window coordinate system conventions.
		// The KHT implementation assumes row-major memory alignment for
		// images. Also, it assumes that the origin of the image coordinate
		// system is at the center of the image, with the x-axis growing to
		// the right and the y-axis growing down.
		if (sin_theta != 0.0) {
			p1.x = -width * 0.5; p1.y = (rho - p1.x * cos_theta) / sin_theta;
			p2.x = width * 0.5 - 1; p2.y = (rho - p2.x * cos_theta) / sin_theta;
		}
		else {
			p1.x = rho; p1.y = -height * 0.5;
			p2.x = rho; p2.y = height * 0.5 - 1;
		}
		p1.x += width * 0.5; p1.y += height * 0.5;
		p2.x += width * 0.5; p2.y += height * 0.5;
		Line(img, p1.x, p1.y, p2.x, p2.y, 0xffffffff);
	}
}

void DrawLinesTrd(STMImage& img, vector <Vec2f> lines)
{
	double rho, theta;
	Point pt1, pt2;
	double a, b;
	for (size_t ii = 0; ii < lines.size(); ii++)
	{
		rho = lines[ii][0]/3.0; theta = lines[ii][1];

		Line(img, rho, theta, 0xffffffff);
	}
}

void DrawLinesTrd(STMImage &img, vector<Vec4i> &lines)
{
	for (size_t j = 0; j < lines.size(); j++)
	{
		Vec4i l = lines[j];
		Line(img, l[0] / 3, l[1] / 3, l[2] / 3, l[3] / 3, 0xffffffff);		
	}
}

void DrawLinesTrd(STMImage& img, kht::ListOfLines lines, int maxlines)
{
	cv::Point p1, p2;
	int width = img.width*3, height = img.height*3;

	for (std::size_t j = 0; j != maxlines; ++j) {
		auto const& line = lines[j];
		std::double_t rho = line.rho;
		std::double_t theta = line.theta * DEGREES_TO_RADIANS;
		std::double_t cos_theta = cos(theta), sin_theta = sin(theta);

		// Convert from KHT to OpenCV window coordinate system conventions.
		// The KHT implementation assumes row-major memory alignment for
		// images. Also, it assumes that the origin of the image coordinate
		// system is at the center of the image, with the x-axis growing to
		// the right and the y-axis growing down.
		if (sin_theta != 0.0) {
			p1.x = -width * 0.5; p1.y = (rho - p1.x * cos_theta) / sin_theta;
			p2.x = width * 0.5 - 1; p2.y = (rho - p2.x * cos_theta) / sin_theta;
		}
		else {
			p1.x = rho; p1.y = -height * 0.5;
			p2.x = rho; p2.y = height * 0.5 - 1;
		}
		p1.x += width * 0.5; p1.y += height * 0.5;
		p2.x += width * 0.5; p2.y += height * 0.5;
		Line(img, p1.x/3, p1.y/3, p2.x/3, p2.y/3, 0xffffffff);
	}
}

int MyCircularHough(Mat &img, vector <Vec3f> &circles)
{
	HoughCircles(img, circles, CV_HOUGH_GRADIENT, 2, img.rows / 6, 300, 100, 5, 60);
	return circles.size();
}

int MyCircularHoughx3(Mat &img, vector <Vec3f> &circles)
{
	HoughCircles(img, circles, CV_HOUGH_GRADIENT, 2, img.rows / 6, 300, 100, 15, 180);
	return circles.size();
}

void DrawCircles(STMImage &img, vector <Vec3f> &circles)
{
	for (size_t ii = 0; ii < circles.size(); ii++)
	{
		Point center(cvRound(circles[ii][0]), cvRound(circles[ii][1]));
		int radius = cvRound(circles[ii][2]);
		Circle(img, center.x, center.y, radius);
		Circle(img, center.x, center.y, radius+1);
		Circle(img, center.x, center.y, radius +2);
	}
}

void DrawCirclesTrd(STMImage &img, vector <Vec3f> &circles)
{
	for (size_t ii = 0; ii < circles.size(); ii++)
	{
		Point center(cvRound(circles[ii][0]), cvRound(circles[ii][1]));
		int radius = cvRound(circles[ii][2]);
		Circle(img, center.x/3, center.y/3, radius/3);
		Circle(img, center.x / 3, center.y / 3, radius / 3+1);
		Circle(img, center.x / 3, center.y / 3, radius / 3+2);
	}
}

void SpecialThin(Mat &edges)
{
	int x, y;
	for (y = 1; y < edges.cols - 1; y++)
		for (x = 1; x < edges.rows - 1; x++)
		{
			if (edges.at<unsigned char>(x, y) > 1)
			{
				if (edges.at<unsigned char>(x - 1, y) > 1 && edges.at<unsigned char>(x, y - 1) > 1) edges.at<unsigned char>(x, y) = 0;
				else if (edges.at<unsigned char>(x - 1, y) > 1 && edges.at<unsigned char>(x, y + 1) > 1) edges.at<unsigned char>(x, y) = 0;
				else if (edges.at<unsigned char>(x + 1, y) > 1 && edges.at<unsigned char>(x, y + 1) > 1) edges.at<unsigned char>(x, y) = 0;
				else if (edges.at<unsigned char>(x + 1, y) > 1 && edges.at<unsigned char>(x, y - 1) > 1) edges.at<unsigned char>(x, y) = 0;
			}
		}
}

void Smooth(Mat &stretch, Mat &tt)
{
	int x, y;
	
	stretch.copyTo(tt);
	uchar *s = stretch.data;
	uchar *t = tt.data;
	int X = tt.cols;
	
	for (x = 1; x < stretch.cols - 1; x++)
		for (y = 1; y < stretch.rows - 1; y++)
		{
			if (stretch.at<unsigned char>(y, x) > 0)
			{				
				//Fig2-1
				//if (stretch.at<unsigned char>(y, x-1) > 0 && stretch.at<unsigned char>(y + 1, x + 1) > 0)
				if(s[X*y+x-1]>0 && s[X*(y+1) + x - 1])
				{
					//tt.at<unsigned char>(x, y) = 0; tt.at<unsigned char>(x, y + 1) = 255;
					t[X*y + x] = 0; t[X*(y + 1) + x] = 255;
				}//Fig2-2
				//else if (stretch.at<unsigned char>(y, x-1) > 0 && stretch.at<unsigned char>(y - 1, x + 1) > 0)
				else if (s[X*y + x - 1] > 0 && s[X*(y -1) + x + 1])
				{
					//tt.at<unsigned char>(x, y) = 0; tt.at<unsigned char>(x, y - 1) = 255;
					t[X*y + x] = 0; t[X*(y - 1) + x] = 255;
				}//Fig2-3
				//else if (stretch.at<unsigned char>(y - 1, x - 1) > 0 && stretch.at<unsigned char>(y, x+1) > 0)
				else if (s[X*(y-1) + x - 1] > 0 && s[X*(y) + x + 1])
				{
					//tt.at<unsigned char>(x, y) = 0; tt.at<unsigned char>(x, y - 1) = 255;
					t[X*y + x] = 0; t[X*(y - 1) + x] = 255;
				}//Fig2-4
				//else if (stretch.at<unsigned char>(y+1, x- 1) > 0 && stretch.at<unsigned char>(y,x + 1) > 0)
				if (s[X*(y+1) + x - 1] > 0 && s[X*(y ) + x + 1])
				{
					//tt.at<unsigned char>(x, y) = 0; tt.at<unsigned char>(x, y + 1) = 255;
					t[X*y + x] = 0; t[X*(y + 1) + x] = 255;
				}//Fig2-5
				//else if (stretch.at<unsigned char>(y - 1, x - 1) > 0 && stretch.at<unsigned char>(y + 1,x) > 0)
				if (s[X*(y-1) + x - 1] > 0 && s[X*(y + 1) + x])
				{
					//tt.at<unsigned char>(x, y) = 0; tt.at<unsigned char>(x - 1, y) = 255;
					t[X*y + x] = 0; t[X*(y ) + x-1] = 255;
				}//Fig2-6
				//else if (stretch.at<unsigned char>(y-1, x+1) > 0 && stretch.at<unsigned char>(y + 1,x) > 0)
				if (s[X*y + x - 1] > 0 && s[X*(y + 1) + x - 1])
				{
					//tt.at<unsigned char>(x, y) = 0; tt.at<unsigned char>(x + 1, y) = 255;
					t[X*y + x] = 0; t[X*(y ) + x+1] = 255;
				}//Fig2-7
				//else if (stretch.at<unsigned char>(y - 1,x) > 0 && stretch.at<unsigned char>(y + 1, x + 1) > 0)
				if (s[X*(y-1) + x] > 0 && s[X*(y + 1) + x + 1])
				{
					//tt.at<unsigned char>(x, y) = 0; tt.at<unsigned char>(x + 1, y) = 255;
					t[X*y + x] = 0; t[X*(y ) + x+1] = 255;
				}//Fig2-8
				//else if (stretch.at<unsigned char>( y - 1,x) > 0 && stretch.at<unsigned char>(y + 1, x - 1) > 0)
				if (s[X*(y-1) + x] > 0 && s[X*(y + 1) + x - 1])
				{
					//tt.at<unsigned char>(x, y) = 0; tt.at<unsigned char>(x - 1, y) = 255;
					t[X*y + x] = 0; t[X*(y ) + x-1] = 255;
				}
			}
		}
	
}
void Stretch(Mat &edges,Mat &stretch)
{
	int x, y;
	stretch = Mat::zeros(edges.rows * 3, edges.cols * 3, edges.type());
	for (y = 1; y < edges.cols - 1; y++)
		for (x = 1; x < edges.rows - 1; x++)
		{
			if (edges.at<unsigned char>(x, y) > 1)
			{
				//ust kenar
				stretch.at<unsigned char>(x * 3 - 3, y * 3 - 3) = edges.at<unsigned char>(x - 1, y - 1);
				stretch.at<unsigned char>(x * 3 - 2, y * 3 - 3) = edges.at<unsigned char>(x, y - 1);
				stretch.at<unsigned char>(x * 3 - 1, y * 3 - 3) = edges.at<unsigned char>(x + 1, y - 1);
				//orta
				stretch.at<unsigned char>(x * 3 - 3, y * 3 - 2) = edges.at<unsigned char>(x - 1, y);
				stretch.at<unsigned char>(x * 3 - 1, y * 3 - 2) = edges.at<unsigned char>(x + 1, y);
				stretch.at<unsigned char>(x * 3 - 2, y * 3 - 2) = edges.at<unsigned char>(x, y);
				//alt
				stretch.at<unsigned char>(x * 3 - 3, y * 3 - 1) = edges.at<unsigned char>(x - 1, y + 1);
				stretch.at<unsigned char>(x * 3 - 2, y * 3 - 1) = edges.at<unsigned char>(x, y + 1);
				stretch.at<unsigned char>(x * 3 - 1, y * 3 + -1) = edges.at<unsigned char>(x + 1, y + 1);
			}
		}
}


LRESULT CALLBACK WndProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	switch (uMsg)
	{

	case WM_CREATE:
	{

		linesbmp = CreateWindow(WC_STATIC, "", WS_CHILD | WS_VISIBLE | SS_BITMAP | WS_THICKFRAME,
			1260, 0, 640, 480, hWnd, NULL, hInst, NULL);
		mesgbox = CreateWindow(WC_EDIT, "", WS_CHILD | WS_VISIBLE | WS_DLGFRAME | WS_VSCROLL | ES_MULTILINE,
			1290, 750, 450, 140, hWnd, NULL, hInst, NULL);

		ProbHoughstatic = CreateWindow(WC_STATIC, "Pr. Probabilistic H. Transform Parameters:", WS_CHILD | WS_VISIBLE,
			1290, 703, 300, 20, hWnd, NULL, hInst, NULL);
		ProbHoughEd1 = CreateWindow(WC_EDIT, "180.0", WS_CHILD | WS_VISIBLE | WS_DLGFRAME,
			1290, 723, 55, 20, hWnd, NULL, hInst, NULL);
		ProbHoughEd2 = CreateWindow(WC_EDIT, "40", WS_CHILD | WS_VISIBLE | WS_DLGFRAME,
			1350, 723, 55, 20, hWnd, NULL, hInst, NULL);
		ProbHoughEd3 = CreateWindow(WC_EDIT, "40.0", WS_CHILD | WS_VISIBLE | WS_DLGFRAME,
			1410, 723, 55, 20, hWnd, NULL, hInst, NULL);
		ProbHoughEd4 = CreateWindow(WC_EDIT, "6.0", WS_CHILD | WS_VISIBLE | WS_DLGFRAME,
			1470, 723, 55, 20, hWnd, NULL, hInst, NULL);

		msg.SetWindow(mesgbox);
		msg.Post("!!! This program requires  HD (1920x1080) resolution monitor!!!\n");
		msg.Post("!!! 4K is even better!!!\n");
	}
	break;
	case WM_COMMAND:
		switch (LOWORD(wParam))
		{
		case IDC_OPEN:
		{
			OPENFILENAME ofn;
			image.zeros(1, 1, CV_8U);
			char szFile[260];
			ZeroMemory(&ofn, sizeof(ofn));
			ofn.lStructSize = sizeof(ofn);
			ofn.hwndOwner = hWnd;
			ofn.lpstrFile = szFile;
			ofn.lpstrFile[0] = '\0';
			ofn.nMaxFile = sizeof(szFile);
			ofn.lpstrFilter = "All\0*.*\0Text\0*.TXT\0";
			ofn.nFilterIndex = 1;
			ofn.lpstrFileTitle = NULL;
			ofn.nMaxFileTitle = 0;
			ofn.lpstrInitialDir = NULL;
			ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;
			GetOpenFileName(&ofn);

			stretch = Mat::zeros(1, 1, CV_8U);
			edges = Mat::zeros(1, 1, CV_8U);
			//Mat del(512, 640, CV_8U);
			//del = 128;
			image = cv::imread(szFile, 0);//if 0 grayscale
			imcolor = cv::imread(szFile, 1);
			if (image.rows < 10) break;
			sprintf_s(mes, 64, "image resolution:%d x %d ", image.rows, image.cols);
			SetWindowText(mesgbox, mes);
			msg.Post("\nThis program is written by Dr. Ibrahim Cem BAYKAL.\n");

			//cv::imshow("My Window", del);
			//cv::imshow("Edge Window", del);
			//ShowWindow(realparent, SW_HIDE);
			UpdateWindow(hWnd);
			MoveWindow(mainbmp, 3, 0, image.cols, image.rows, TRUE);
			cv::imshow("My Window", image);
			ShowWindow(realparent, SW_HIDE);
			copyofimg = image;
			blinearflag = 0;

			//-------------------------------------------------------------
			HDC lineDC = GetDC(linesbmp);
			i.width = image.cols;
			i.height = image.rows;
			FillBMPInf(&i.info, i.width, i.height, 32);
			unsigned char u;
			i.HBMRGB = CreateDIBSection(NULL, &i.info, DIB_RGB_COLORS, (void **)&i.RGB, 0, 0);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
			//--------------------------------
			
		}
		break;
		//ORIGINAL METHOD --------------------------------------------------------------------
		case IDC_CANNY:
		{
			if (image.rows < 3)break;			
			Mat temp;
			GaussianBlur(image, temp, Size(7, 7), 2.1, 2.1);
			Canny(temp, edges, 18, 40, 3);
			MoveWindow(mainbmp, 3, 0, edges.cols, edges.rows, TRUE);
			cv::imshow("My Window", edges);
		}
		break;
		case IDC_STD_HOUGH:
		{
			vector <Vec2f> lines;
			int numberoflines = MyStdHough(edges, lines);
			DrawLines(i, lines);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
			itoa(numberoflines, tmes, 10);
			//strcpy_s(mes, 64, "Number of Lines: ");
			//strcat_s(mes, 64, tmes);
			//SetWindowText(mesgbox, mes);
		}
		break;
		case IDC_STD_KHT:
		{
			cv::Point p1, p2;
			kht::ListOfLines lines;
			Mat tmp;
			int width= edges.cols, height= edges.rows;
			tmp = edges;
			kht::run_kht(lines, tmp.ptr(), width, height);
			CopyToSTM(copyofimg, i);
			DrawLines(i, lines,30);			
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);			
		}
		break;
		case IDC_PRO_HOUGH:
		{
			vector<Vec4i> lines;
			int numberoflines=MyHough(edges,lines);
			DrawLines(i, lines);			
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
			itoa(numberoflines, tmes, 10);
			strcpy_s(mes, 64, "Number of Lines: ");
			strcat_s(mes, 64, tmes);
			SetWindowText(mesgbox, mes);
		}
		break;
		case IDC_CRCL_HOUGH:
		{
			vector <Vec3f> circles;
			MyCircularHough(edges, circles);
			CopyToSTM(copyofimg, i);
			//ColorCopy(imcolor, i);
			DrawCircles(i, circles);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
		}
		break;
		//BILINEAR TRANSFORM -------------------------------------------------------
		case IDC_BILINEAR:
		{
			cv::resize(image, image_x3, Size(), 3.0, 3.0, INTER_LINEAR);
			MoveWindow(mainbmp, 3, 0, image_x3.cols, image_x3.rows, TRUE);
			cv::imshow("My Window", image_x3);
			blinearflag = 1;
		}
		break;
		case IDC_CANNY_x3:
		{
			if (image.rows < 3)break;
			stretch=Mat::zeros(1, 1, CV_8U);
			edges=Mat::zeros(1, 1, CV_8U);
			Mat temp;
			GaussianBlur(image_x3, temp, Size(13, 13), 2.1, 2.1);
			Canny(temp, edges_x3, 18, 40, 3);
			MoveWindow(mainbmp, 3, 0, edges_x3.cols, edges_x3.rows, TRUE);
			cv::imshow("My Window", edges_x3);
		}
		break;
		case IDC_STD_HOUGH_BIL:
		{
			vector <Vec2f> lines;			
			int numberoflines = MyStdHoughx3(edges_x3, lines);
			DrawLinesTrd(i, lines);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
			itoa(numberoflines, tmes, 10);
			strcpy_s(mes, 64, "Number of Lines: ");
			strcat_s(mes, 64, tmes);
			SetWindowText(mesgbox, mes);
		}
		break;
		case IDC_KHT_BIL:
		{
			cv::Point p1, p2;
			kht::ListOfLines lines;
			Mat tmp;
			int width = edges_x3.cols, height = edges_x3.rows;
			tmp = edges_x3;
			kht::run_kht(lines, tmp.ptr(), width, height);			
			CopyToSTM(copyofimg, i);
			DrawLinesTrd(i, lines, 30);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
		}
		break;
		case IDC_PRO_HOUGH_BIL:
		{
			vector<Vec4i> lines;
			int numberoflines = MyHoughx3(edges_x3, lines);
			DrawLinesTrd(i, lines);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
			itoa(numberoflines, tmes, 10);
			strcpy_s(mes, 64, "Number of Lines: ");
			strcat_s(mes, 64, tmes);
			SetWindowText(mesgbox, mes);
		}
		break;
		case IDC_CRCL_HOUGH_BIL:
		{
			vector <Vec3f> circles;
			MyCircularHoughx3(edges_x3, circles);
			CopyToSTM(copyofimg, i);
			//ColorCopy(imcolor, i);
			DrawCirclesTrd(i, circles);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
			Sleep(10);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
		}
		break;
		// MORPHOLOGIC STRETCH -----------------------------------------------------------------
		case IDC_THIN:
		{
			if (edges.rows < 3)break;
			stretch = Mat::zeros(1, 1, CV_8U);
			SpecialThin(edges);
			
			MoveWindow(mainbmp, 3, 0, edges.cols, edges.rows, TRUE);
			cv::imshow("My Window", edges);
		}
		break;
		case IDC_STRETCH:
		{
			if (edges.rows < 3)break;
			Stretch(edges, stretch);
			MoveWindow(mainbmp, 3, 0, stretch.cols, stretch.rows, TRUE);
			cv::imshow("My Window", stretch);
		}
		break;
		case IDC_SMOOTH:
		{	
			if (stretch.rows < 3)break;
			Smooth(stretch,smooth);
			MoveWindow(mainbmp, 3, 0, smooth.cols, smooth.rows, TRUE);
			cv::imshow("My Window", smooth);
		}
		break;
		case IDC_STD_HOUGH_STRTCH:
		{
			vector <Vec2f> lines;
			int numberoflines = MyStdHoughx3(smooth, lines);
			DrawLinesTrd(i, lines);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
			itoa(numberoflines, tmes, 10);
			strcpy_s(mes, 64, "Number of Lines: ");
			strcat_s(mes, 64, tmes);
			SetWindowText(mesgbox, mes);
		}
		break;
		case IDC_KHT_STRTCH:
		{
			cv::Point p1, p2;
			kht::ListOfLines lines;
			Mat tmp;
			int width = smooth.cols, height = smooth.rows;
			tmp = smooth;	
			kht::run_kht(lines, tmp.ptr(), width, height);
			CopyToSTM(copyofimg, i);
			DrawLinesTrd(i, lines, 30);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
		}
		break;
		case IDC_PRO_HOUGH_STRTCH:
		{
			vector<Vec4i> lines;
			int numberoflines = MyHoughx3(smooth, lines);
			DrawLinesTrd(i, lines);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
			itoa(numberoflines, tmes, 10);
			strcpy_s(mes, 64, "Number of Lines: ");
			strcat_s(mes, 64, tmes);
			SetWindowText(mesgbox, mes);
		}
		break;		
		case IDC_CRCL_HOUGH_STRTCH:
		{
			vector <Vec3f> circles;
			MyCircularHoughx3(smooth, circles);
			CopyToSTM(copyofimg, i);
			//ColorCopy(imcolor, i);
			DrawCirclesTrd(i, circles);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
		}
		break;
		
		//---------------------------------------- TESTS -----------------------------
		case IDC_BILIN_STD_HOUGH_TEST:
		{
			Mat temp;
			vector <Vec2f> lines;
			
			double duration;
			chrono::duration<double> elapsed;
			msg.Post("Standard Hough using Bilinear Stretch operator benchmark started...\n");

			auto difstart = chrono::high_resolution_clock::now();
			for (int x = 0; x < 100; x++)
			{
				cv::resize(image, image_x3, Size(), 3.0, 3.0, INTER_LINEAR);
				GaussianBlur(image_x3, temp, Size(13, 13), 2.1, 2.1);
				Canny(temp, edges_x3, 18, 40, 3);
				int numberoflines = MyStdHoughx3(edges_x3, lines);
			}
			auto difend = chrono::high_resolution_clock::now();

			DrawLinesTrd(i, lines);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
			msg.Post("finished.\n");
			elapsed = difend - difstart;
			duration = elapsed.count();
			msg.Post("%f sec. elapsed for 100 iterations\n", duration);
		}
		break;
		case IDC_STRETCH_STD_HOUGH_TEST:
		{
			Mat temp;
			vector <Vec2f> lines;
			double duration;
			chrono::duration<double> elapsed;
			msg.Post("Standard Hough using Morphological stretch operator benchmark started...\n");
			auto difstart = chrono::high_resolution_clock::now();
			for (int x = 0; x < 100; x++)
			{
				GaussianBlur(image, temp, Size(7, 7), 2.1, 2.1);
				Canny(temp, edges, 18, 40, 3);
				SpecialThin(edges);
				Stretch(edges, stretch);
				Smooth(stretch, smooth);
				MyStdHoughx3(smooth, lines);
			}
			auto difend = chrono::high_resolution_clock::now();
			DrawLinesTrd(i, lines);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
			msg.Post("finished.\n");
			elapsed = difend - difstart;
			duration = elapsed.count();
			msg.Post("%f sec. elapsed for 100 iterations\n", duration);

		}
		break;
		case IDC_BILIN_KHT_TEST:
		{
			Mat temp;
			cv::Point p1, p2;
			kht::ListOfLines lines;
			double duration;
			chrono::duration<double> elapsed;
			msg.Post("Kernel Based Hough using Bilinear Stretch operator benchmark started...\n");
			
			
						
			stretch = Mat::zeros(1, 1, CV_8U);
			edges = Mat::zeros(1, 1, CV_8U);

			auto difstart = chrono::high_resolution_clock::now();
			for (int x = 0; x < 100; x++)
			{
				cv::resize(image, image_x3, Size(), 3.0, 3.0, INTER_LINEAR);
				GaussianBlur(image_x3, temp, Size(13, 13), 2.1, 2.1);
				Canny(temp, edges_x3, 18, 40, 3);
				
				int width = edges_x3.cols, height = edges_x3.rows;
				
				kht::run_kht(lines, edges_x3.ptr(), width, height);
			}
			auto difend = chrono::high_resolution_clock::now();
			msg.Post("finished.\n");
			elapsed = difend - difstart;
			duration = elapsed.count();
			msg.Post("%f sec. elapsed for 100 iterations\n", duration);

			CopyToSTM(copyofimg, i);
			DrawLinesTrd(i, lines, 30);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
		}
		break;
		case IDC_STRETCH_KHT_TEST:
		{
			Mat temp;
			cv::Point p1, p2;
			kht::ListOfLines lines;
			double duration;
			chrono::duration<double> elapsed;
			msg.Post("Kernel Based Hough using Bilinear Stretch operator benchmark started...\n");
			
			
			
			stretch = Mat::zeros(1, 1, CV_8U);
			edges = Mat::zeros(1, 1, CV_8U);

			auto difstart = chrono::high_resolution_clock::now();
			for (int x = 0; x < 100; x++)
			{

				GaussianBlur(image, temp, Size(7, 7), 2.1, 2.1);
				Canny(temp, edges, 18, 40, 3);
				SpecialThin(edges);
				Stretch(edges, stretch);
				Smooth(stretch, smooth);
				int width = smooth.cols, height = smooth.rows;

				kht::run_kht(lines, smooth.ptr(), width, height);
			}
			auto difend = chrono::high_resolution_clock::now();
			msg.Post("finished.\n");
			elapsed = difend - difstart;
			duration = elapsed.count();
			msg.Post("%f sec. elapsed for 100 iterations\n", duration);

			CopyToSTM(copyofimg, i);
			DrawLinesTrd(i, lines, 30);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
		}
		break;
		case IDC_BILIN_PROB_HOUGH_TEST:
		{
			Mat temp;
			vector<Vec4i> lines;
			double duration;
			chrono::duration<double> elapsed;
			msg.Post("Prob. Hough using Bilinear Stretch operator benchmark started...\n");
			auto difstart = chrono::high_resolution_clock::now();
			for (int x = 0; x < 100; x++)
			{
				cv::resize(image, image_x3, Size(), 3.0, 3.0, INTER_LINEAR);
				GaussianBlur(image_x3, temp, Size(13, 13), 2.1, 2.1);
				Canny(temp, edges_x3, 18, 40, 3);
				int numberoflines = MyHoughx3(edges_x3, lines);
			}
			auto difend = chrono::high_resolution_clock::now();
			DrawLinesTrd(i, lines);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
			msg.Post("finished.\n");
			elapsed = difend - difstart;
			duration = elapsed.count();
			msg.Post("%f sec. elapsed for 100 iterations\n", duration);
		}
		break;
		case IDC_STRETCH_PROB_HOUGH_TEST:
		{
			Mat temp;
			vector<Vec4i> lines;
			double duration;
			chrono::duration<double> elapsed;
			msg.Post("Probabilistic Hough using orphological stretch operator benchmark started...\n");
			auto difstart = chrono::high_resolution_clock::now();
			for (int x = 0; x < 100; x++)
			{				
				GaussianBlur(image, temp, Size(7, 7), 2.1, 2.1);
				Canny(temp, edges, 18, 40, 3);
				SpecialThin(edges);
				Stretch(edges, stretch);
				Smooth(stretch, smooth);
				MyHoughx3(smooth, lines);
			}
			auto difend = chrono::high_resolution_clock::now();
			DrawLinesTrd(i, lines);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
			msg.Post("finished.\n");
			elapsed = difend - difstart;
			duration = elapsed.count();
			msg.Post("%f sec. elapsed for 100 iterations\n", duration);

		}
		break;
		case IDC_BILIN_CIRC_HOUGH_TEST:
		{
			Mat temp;
			vector <Vec3f> circles;
			double duration;
			chrono::duration<double> elapsed;
			msg.Post("Circular Hough using Bilinear Trsfm benchmark started...\n");
			auto difstart = chrono::high_resolution_clock::now();
			for (int x = 0; x < 10; x++)
			{
				cv::resize(image, image_x3, Size(), 3.0, 3.0, INTER_LINEAR);
				GaussianBlur(image_x3, temp, Size(13, 13), 2.1, 2.1);
				Canny(temp, edges_x3, 18, 40, 3);

				MyCircularHoughx3(edges_x3, circles);

			}
			auto difend = chrono::high_resolution_clock::now();
			CopyToSTM(copyofimg, i);
			DrawCirclesTrd(i, circles);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
			msg.Post("finished.\n");
			elapsed = difend - difstart;
			duration = elapsed.count();
			msg.Post("%f sec. elapsed for 10 iterations\n", duration);
		}
		break;
		case IDC_STRETCH_CIRC_HOUGH_TEST:
		{
			Mat temp;
			vector <Vec3f> circles;
			double duration;
			chrono::duration<double> elapsed;
			msg.Post("Circular Hough Trnsfm using morphological stretch operator benchmark started...\n");
			auto difstart = chrono::high_resolution_clock::now();
			for (int x = 0; x < 10; x++)
			{
				GaussianBlur(image, temp, Size(7, 7), 2.1, 2.1);
				Canny(temp, edges, 18, 40, 3);
				SpecialThin(edges);
				Stretch(edges, stretch);
				Smooth(stretch, smooth);
				MyCircularHoughx3(smooth, circles);

			}
			auto difend = chrono::high_resolution_clock::now();
			CopyToSTM(copyofimg, i);
			DrawCirclesTrd(i, circles);
			SendMessage(linesbmp, STM_SETIMAGE, 0, (LPARAM)i.HBMRGB);
			msg.Post("finished.\n");
			elapsed = difend - difstart;
			duration = elapsed.count();
			msg.Post("%f sec. elapsed for 10 iterations\n", duration);
		}
		break;
		case IDC_HELP:
		{
			MessageBox(NULL, "!!! This program requires  HD (1920x1080) resolution!!!\n 4K is even better.\nFrom the Menu select FILE->Open,\
				\nChoose a menu item \nThen Canny Edge\nPerform algorithms in order from top","HOW TO USE:", MB_OK);
		}
		break;
		case IDC_ABOUT:
		{
			MessageBox(NULL, "This program is written by:\nDr. Ibrahim Cem BAYKAL\nicbaykal@atu.edu.tr", "AUTHOR:", MB_OK);
		}
		break;
		default:
			break;
		}
		break;

	case WM_DESTROY:

		i.~STMImage();
		PostQuitMessage(0);
		break;

	default:
		return(DefWindowProc(hWnd, uMsg, wParam, lParam));
	}

	return(0L);
}
