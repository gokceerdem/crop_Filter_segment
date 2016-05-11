// TestDbProject.cpp : Defines the entry point for the console application.
//
#include <stdio.h>
#include "sqlite3.h"
#include "stdafx.h"
#include <cstdlib>
#include "HEADER.h"

using namespace std;
using namespace cv;

int labelColors[10000][2] = { 0 };
int myClrs[TOTAL_COLOR_Nr + 1][3];
int dummy[10000][3];

vector<Point> contourArray;
vector<vector<Point>> imageContours;
vector<int> segmentClrs;

list<hsvc> col_hash_map;
coor c;

list<SegmentFeatures> featureVector;
SegmentFeatures mF;
list<SegmentMoments> momentVector;
SegmentMoments SegMom;

ImData myImData;
time_t tstart, tend, mstart, mend;


stringstream imgstream, mystream, seq_stream, imgstr, segmentstream, filterstream, cropstream, filtercropstream;

//Static array
int label[pyrHeight][pyrWidth] = { 0 };
int I[pyrHeight][pyrWidth] = { 0 };
int Q[pyrHeight][pyrWidth + 1] = { 0 };

int EQ[MAX_EXP_NrOf_LABELS] = { 0 };
int relationLUT[TOTAL_COLOR_Nr + 1][TOTAL_COLOR_Nr + 1] = { 0 };
int labelNr;

int centerSum[2][1] = { 0 };


static int callback(void *NotUsed, int argc, char **argv, char **azColName){
	int i;
	for (i = 0; i<argc; i++){
		printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
	}
	printf("\n");
	return 0;
}

void LUT(int relationLUT[TOTAL_COLOR_Nr + 1][TOTAL_COLOR_Nr + 1]){


	for (int t = 0; t <TOTAL_COLOR_Nr + 1; t++){
		relationLUT[t][t] = 1;
	}
	relationLUT[24][25] = 1;//mx white & white
	relationLUT[25][24] = 1;//white & mx white 
	//relationLUT[23][22]=1;//black & mx black
	//relationLUT[22][23]=1;
	relationLUT[23][18] = 1;//black & dark gray
	relationLUT[18][23] = 1;
	relationLUT[18][22] = 1;//dark gray & mx black
	relationLUT[22][18] = 1;
	//relationLUT[24][19]=1;//mx white & light gray 
	relationLUT[8][9] = 1; // light blue & dark blue  
	relationLUT[21][3] = 1;//yellow & brown
	relationLUT[3][20] = 1;

}

void colorImg(vector<Mat> channels, Mat ch0, Mat ch1, Mat ch2, int I[pyrHeight][pyrWidth], int imgnr){
	channels.clear();

	myColors(myClrs);

	for (int colorindex = 1; colorindex<TOTAL_COLOR_Nr + 1; colorindex++)
	{
		for (int i = 0; i < myImData.h; i++) {
			for (int j = 0; j < myImData.w; j++) {

				int colorval = I[i][j];

				if (colorval == colorindex){
					ch2.at<uchar>(i, j) = myClrs[colorindex][0];
					ch1.at<uchar>(i, j) = myClrs[colorindex][1];
					ch0.at<uchar>(i, j) = myClrs[colorindex][2];
				}

			}
		}
	}
	Mat fin_img2;
	channels.push_back(ch0); //b
	channels.push_back(ch1); //g
	channels.push_back(ch2); //r
	merge(channels, fin_img2);

	//imgstr << "colored_by_respected_colors_"<< imagename << quant << type;
	imgstr << imagename1 << "_" << imgnr << quant << type;
	string myfilename = imgstr.str();
	imgstr.str("");
	imwrite(myfilename, fin_img2);

	ch0.release();
	ch1.release();
	ch2.release();
	fin_img2.release();
	channels.clear();
}

void Labeling(int &labelNr, int label[pyrHeight][pyrWidth], int I[pyrHeight][pyrWidth], int Q[pyrHeight][pyrWidth + 1], int EQ[MAX_EXP_NrOf_LABELS], ImData &myImData){
	// Label (0,0) start point
	int L = 0;
	labelNr = 0;
	++L; ++labelNr;
	EQ[L] = (L);
	label[0][0] = L; Q[0][1] = L;

	// Label first row 	
	for (int x = 1; x<myImData.w; x++){

		int y = 0;
		int n1x = x - 1;

		if (I[y][n1x] == I[y][x]){
			label[y][x] = label[y][n1x];
			Q[y][x + 1] = label[y][x];
		}
		if (I[y][n1x] != I[y][x]){
			++L; ++labelNr;
			EQ[L] = (L);
			label[y][x] = L;
			Q[y][x + 1] = L;
		}
	}

	// Label first column starting from second row	

	for (int y = 1; y<myImData.h; y++){
		if (I[y][0] == I[y - 1][0]){
			label[y][0] = label[y - 1][0];
			Q[y][1] = label[y][0];
		}

		if (I[y][0] != I[y - 1][0]){
			++L; ++labelNr;
			EQ[L] = (L);
			label[y][0] = L;
			Q[y][1] = label[y][0];
		}
	}

	//Label the rest of the img

	for (int x = 1; x<myImData.w; x++){
		for (int y = 1; y<myImData.h; y++){


			int sx = x - 1; int sy = y;
			int tx = x;	 int ty = y - 1;

			if (I[y][x] == I[sy][sx] && I[y][x] != I[ty][tx]){
				label[y][x] = label[sy][sx];
			}

			if (I[y][x] != I[sy][sx] && I[y][x] == I[ty][tx]){
				label[y][x] = label[ty][tx];
			}

			if (I[y][x] != I[sy][sx] && I[y][x] != I[ty][tx]){
				++L; ++labelNr;
				EQ[L] = (L);
				label[y][x] = L;
			}

			if (I[y][x] == I[sy][sx] && I[y][x] == I[ty][tx] && label[sy][sx] == label[ty][tx]){
				label[y][x] = label[ty][tx];
			}

			if (I[y][x] == I[sy][sx] && I[y][x] == I[ty][tx] && label[sy][sx] != label[ty][tx]){
				int comp = (label[sy][sx]<label[ty][tx]); // Ls < Lt -->  1
				// Ls > Lt -->  0
				int L1, L2; //L1<L2
				comp ? L1 = label[sy][sx] : L1 = label[ty][tx];
				comp ? L2 = label[ty][tx] : L2 = label[sy][sx];

				label[y][x] = L1;
				EQ[L2] = L1;
			}
			Q[y][x + 1] = label[y][x];

		}
	}
	for (int i = 0; i<myImData.h; i++){
		Q[i][0] = label[i][1];

	}

}

void LabelEqualization(int EQ[MAX_EXP_NrOf_LABELS], int label[pyrHeight][pyrWidth], ImData &myImData, int labelColors[10000][2]){
	//Equalization of labels
	for (int k = 1; k<MAX_EXP_NrOf_LABELS; k++){

		if (EQ[k] == 0){ break; }

		if (EQ[k] != k){
			EQ[EQ[k]] == EQ[k] ? 1 : EQ[k] = EQ[EQ[EQ[k]]];
		}

		for (int i = 0; i < myImData.h; i++) {
			for (int j = 0; j < myImData.w; j++) {
				if (label[i][j] == k)
					label[i][j] = EQ[k];
				Q[i][j + 1] = label[i][j];

				//labelColors[label[i][j]][0] = I[i][j];
			}
		}
	}
}

void createHash(String dy){

	string line;
	ifstream myfile(dy);

	while (getline(myfile, line))
	{
		stringstream   linestream(line);
		string         data;
		int hl, hh, sl, sh, vl, vh, color_name;// HueLow, HueHigh, SaturationLow, SaturationHigh, ValueLow, ValueHigh

		getline(linestream, data, '\t');

		linestream >> hh >> sl >> sh >> vl >> vh >> color_name;
		hsvc new_hsvc;
		new_hsvc.hlow = atoi(data.c_str());
		new_hsvc.hhigh = hh;
		new_hsvc.slow = sl;
		new_hsvc.shigh = sh;
		new_hsvc.vlow = vl;
		new_hsvc.vhigh = vh;
		new_hsvc.col_name = color_name;

		col_hash_map.push_back(new_hsvc);

	}
	myfile.close();

}

void keepcolors(String clrs){
	int cidx = 0;
	string line;
	ifstream file(clrs);

	while (getline(file, line))
	{
		stringstream   linestream(line);
		string        data;
		int r, g, b;

		getline(linestream, data, '\t');
		linestream >> g >> b;

		/* Array implementation*/
		dummy[cidx][0] = atoi(line.c_str());
		dummy[cidx][1] = g;
		dummy[cidx][2] = b;

		++cidx;

	}
	file.close();

}

Mat readFilter(int filtreSirasi, int filterSize)
{

	Mat filterOrg = Mat(filterSize, filterSize, CV_32FC1);
	Mat filter;

	// QFile file(dirr);
	char filtreDosya[50];
	sprintf_s(filtreDosya, "filtreler\\filtre%i.txt", filtreSirasi);
	fstream file_op(filtreDosya, ios::in);

	for (int j = 0; j<filterSize; j++)
		for (int k = 0; k<filterSize; k++){
			file_op >> filterOrg.at<float>(k, j);
		}

	file_op.close();

	convertScaleAbs(filterOrg, filter, 128, 128);

	Mat resizedFilter;

	resize(filter, resizedFilter, resizedFilter.size(), 5, 5);

	//namedWindow("filter");

	//cv::imshow("filter",resizedFilter);

	//waitKey();

	//destroyWindow("filter");


	return filterOrg;
}

void applyFilter(ImData &mid, Mat filter_out, int index) //this function gets the address of ImData and routes filters directly into filtered outputs in this struct.
{
	Mat result = Mat::zeros(mid.h, mid.w, CV_8UC1);
	mid.filter.push_back(Mat::zeros(mid.h, mid.w, CV_8UC1));

	GaussianBlur(mid.intensity, mid.filter[index], Size(5, 5), 5, 5);

	filter2D(mid.filter.at(index), result, result.depth(), filter_out);

	mid.filter[index] = result;
}


void hesapla(int, int, int);



int main()
{
	LUT(relationLUT);
	keepcolors("SegmentColors.txt");
	createHash("ColorQuantas.txt");

	char choice;
	do{
		int image_no, segment_no, filter_no;
		cout << "Image no (10-106): ";
		cin >> image_no;
		cout << "Segment no : (max 150)";
		cin >> segment_no;
		cout << "Filter no : (max 54)";
		cin >> filter_no;



		hesapla(image_no, segment_no, filter_no);
		cout << "Enter choice : (C-Continiue//E-Exit)" << endl;
		cin >> choice;

	} while (choice == 'C' || choice == 'c');



	return 0;
}


void hesapla(int image_no, int segment_no, int filter_no){

	int ivar = image_no;


	mystream << imagename1 << std::setfill('0') << std::setw(3) << ivar << type;
	string myfilename = mystream.str();
	mystream.str("");
	Mat src = imread(myfilename, CV_LOAD_IMAGE_COLOR);

	Mat dst;
	pyrDown(src, dst, Size(src.cols / 2, src.rows / 2));
	bilateralFilter(dst, myImData.original, 9, 30, 30);
	//	bilateralFilter (src,myImData.original, 9, 30, 30 );

	vector<Mat> dst_channels;
	split(myImData.original, dst_channels);

	cvtColor(myImData.original, myImData.intensity, CV_BGR2GRAY);
	cvtColor(myImData.original, myImData.hsvImg, CV_BGR2HSV);

	myImData.h = myImData.original.rows;
	myImData.w = myImData.original.cols;

	vector<Mat> hsvchannels;
	split(myImData.hsvImg, hsvchannels);

	myImData.hsv_filter.push_back(hsvchannels[0]);
	myImData.hsv_filter.push_back(hsvchannels[1]);
	myImData.hsv_filter.push_back(hsvchannels[2]);


	hsvchannels.clear();



	//keep filter outputs in a new subfolder under current path
	// "FilterOutputs" is the name of the folder that will contain filter outputs of the images
	// imagename_filter+filternumber is the general organization of the output images
	string myfolderName = "FilterOutputs";
	string folderCreateCommand = "mkdir " + myfolderName;
	system(folderCreateCommand.c_str());

	//keep filters in here
	vector<Mat> filters;

	for (int f = 0; f<54; f++)
	{
		filters.push_back(readFilter(f, 29));
	}
	for (int i = 0; i<54; i++){
		//filter original image and create filtered outputs
		applyFilter(myImData, filters[i], i);

		/*	//visualise filter outputs
		//imshow("filter",myImData.filter.at(i));
		//waitKey(0);
		//destroyWindow("filter");
		*/

		/* //Save filter outputs in current path
		filterstream << imagename1 <<std::setfill('0') << std::setw(3) << ivar << "filter"<<i<< type;
		string filteroutput = filterstream.str();
		filterstream.str("");
		imwrite(filteroutput,myImData.filter.at(i));
		*/

		filterstream << myfolderName << "/" << imagename1 << std::setfill('0') << std::setw(3) << ivar << "_filter" << i << type;
		string filteroutput = filterstream.str();
		filterstream.str("");
		imwrite(filteroutput, myImData.filter.at(i));
	}
	myImData.filter.push_back(myImData.hsv_filter.at(0)); //Hue filter result added to 5 spatial filter result
	filters.clear();
	//int count[TOTAL_COLOR_Nr];
	//std::fill( count, count+TOTAL_COLOR_Nr, 0);


	for (int k = 0; k < (myImData.h * myImData.w); k++){
		int x = k % (myImData.w);
		int y = (k - x) % (myImData.w - 1);
		bool flag = false;

		for (list<hsvc>::iterator it = col_hash_map.begin(); it != col_hash_map.end(); ++it){


			int val_h = myImData.hsv_filter.at(0).at<uchar>(y, x);
			int val_s = myImData.hsv_filter.at(1).at<uchar>(y, x);
			int val_v = myImData.hsv_filter.at(2).at<uchar>(y, x);


			if (val_h >= it->hlow && val_h <= it->hhigh && val_s >= it->slow &&
				val_s <= it->shigh && val_v >= it->vlow && val_v <= it->vhigh){

				I[y][x] = it->col_name;

				flag = true; break;
			}

		}
	}

	Labeling(labelNr, label, I, Q, EQ, myImData);

	LabelEqualization(EQ, label, myImData, labelColors);


	/*Merge small components with their nearest component*/

	std::unordered_map<int, int> occurrences;

	for (int i = 0; i < myImData.h; ++i){
		for (int j = 0; j < myImData.w; ++j){

			++occurrences[label[i][j]];
		}
	}

	for (int i = 0; i < myImData.h; ++i){
		for (int j = 0; j < myImData.w; ++j){

			if (occurrences[label[i][j]] < MAX_PxNr_SMALL_AREA) {
				EQ[label[i][j]] = Q[i][j];
				Q[i][j + 1] = Q[i][j];

			}
		}
	}
	occurrences.clear();

	// LabelEqualization(EQ, label, myImData,labelColors);

	vector<int> nIndx;
	int indx = 1;

	while (indx != labelNr + 1){
		contourArray.clear();
		for (int i = 0; i < myImData.h; i++) {
			for (int j = 0; j < myImData.w; j++) {

				int val = label[i][j];
				if (val == indx){
					contourArray.push_back(Point(j, i));
				}
			}
		}
		if (contourArray.empty() == false){

			myImData.connComp.push_back(contourArray);
			int clrv;

			// int clrv = I[contourArray.at(0).y][contourArray.at(0).x]; //yanlýþ renk birleþtirme olmasýn diye deðiþtiriyorum burayý 07.10.2015
			if (contourArray.size() == 1){
				clrv = I[contourArray.at(0).y][contourArray.at(0).x];
			}
			if (contourArray.size()>1){
				clrv = I[contourArray.at(contourArray.size() - 1).y][contourArray.at(contourArray.size() - 1).x];
			}

			segmentClrs.push_back(clrv);

			if (contourArray.size()> MIN_PxNr_BIG_AREA){
				nIndx.push_back(myImData.connComp.size() - 1);
			}

		}
		++indx;
	}

	for (int nfc = 0; nfc< nIndx.size(); nfc++){

		int numberofcomponents = nIndx.at(nfc);
		Mat component_Img = Mat::zeros(myImData.h, myImData.w, CV_8UC1);
		Mat dilated_component_Img, dst;
		Mat eroded_;
		// Create binary image of big segment
		for (int comp = 0; comp < myImData.connComp.at(numberofcomponents).size(); comp++){
			Point component = myImData.connComp.at(numberofcomponents).at(comp);
			component_Img.at<uchar>(component.y, component.x) = 255;
		}

		dilate(component_Img, dilated_component_Img, dilation_element);

		// Obtain adjacent parts
		cv::bitwise_xor(component_Img, dilated_component_Img, dst);
		//imshow("dst", dst); //adjacent img
		//	waitKey(0);


		vector<Point> nonZeroCoordinates;		//keep adjacent pixels in here
		findNonZero(dst, nonZeroCoordinates);

		int	ColorNr1 = segmentClrs.at(numberofcomponents);
		int ColorNr2;
		int newLabel = label[myImData.connComp.at(numberofcomponents).at(0).y][myImData.connComp.at(numberofcomponents).at(0).x];

		for (int g = 0; g<nonZeroCoordinates.size(); g++){

			Point AdjPoint = nonZeroCoordinates.at(g);
			ColorNr2 = I[AdjPoint.y][AdjPoint.x];
			if (relationLUT[ColorNr1][ColorNr2] == 1){

				EQ[label[AdjPoint.y][AdjPoint.x]] = newLabel;

			}

		}

	}

	LabelEqualization(EQ, label, myImData, labelColors);


	/*	To see the most current components again push back components and visualise */
	vector<Mat> mychannels;
	Mat fin_img;

	//Mat ch0 =  Mat::zeros(myImData.h,myImData.w,CV_8UC1); //B
	//Mat ch1 =  Mat::zeros(myImData.h,myImData.w,CV_8UC1); //G 
	//Mat ch2 =  Mat::zeros(myImData.h,myImData.w,CV_8UC1); //R

	myImData.connComp.clear();

	int nindx = 1;
	while (nindx != labelNr + 1){

		contourArray.clear();
		for (int i = 0; i < myImData.h; i++) {
			for (int j = 0; j < myImData.w; j++) {

				int nval = label[i][j];
				if (nval == nindx){

					contourArray.push_back(Point(j, i));

					//		ch2.at<uchar>(i,j)= dummy[(nindx)][0];
					//		ch1.at<uchar>(i,j)= dummy[(nindx)][1];
					//		ch0.at<uchar>(i,j)= dummy[(nindx)][2];

				}
			}
		}
		if (contourArray.empty() == false){
			myImData.connComp.push_back(contourArray);
		}
		++nindx;
	}



	//	for(int ms= 0; ms<myImData.connComp.size();ms++){
	//		Mat c_Img=Mat::zeros(myImData.h,myImData.w,CV_8UC1);
	//	for (int d =0; d<myImData.connComp.at(ms).size(); d++){
	//		Point component = myImData.connComp.at(ms).at(d);
	//		c_Img.at<uchar>(component.y,component.x)=255;
	//
	//	}
	//			imshow("c_Img", c_Img); 
	//				
	//			cout<<I[myImData.connComp.at(ms).at(0).y][myImData.connComp.at(ms).at(0).x]<<endl;	waitKey(0);
	//}



	vector<Mat> showsegment;
	Mat segmentImage;
	Mat goster_r(myImData.h, myImData.w, CV_8UC1, Scalar(background));
	Mat goster_g(myImData.h, myImData.w, CV_8UC1, Scalar(background));
	Mat goster_b(myImData.h, myImData.w, CV_8UC1, Scalar(background));

	Mat goster_filtre(myImData.h, myImData.w, CV_8UC1, Scalar(background));


	for (int d = 0; d<myImData.connComp.at(segment_no).size(); d++){
		Point component = myImData.connComp.at(segment_no).at(d);


		goster_r.at<uchar>(component.y, component.x) = dst_channels[0].at<uchar>(component.y, component.x);
		goster_g.at<uchar>(component.y, component.x) = dst_channels[1].at<uchar>(component.y, component.x);
		goster_b.at<uchar>(component.y, component.x) = dst_channels[2].at<uchar>(component.y, component.x);


		goster_filtre.at<uchar>(component.y, component.x) = myImData.filter.at(filter_no).at<uchar>(component.y, component.x);
	}

	showsegment.push_back(goster_r);
	showsegment.push_back(goster_g);
	showsegment.push_back(goster_b);

	merge(showsegment, segmentImage);

	imshow("Required Segment", segmentImage);
	waitKey(0);

	destroyWindow("Required Segment");

	imshow("Required Filter Result", goster_filtre);
	waitKey(0);
	destroyWindow("Required Filter Result");


	cropstream << crop << ivar << myzero << segment_no << type;
	string segmentfilename = cropstream.str();
	cropstream.str("");
	imwrite(segmentfilename, segmentImage);

	string CmyfolderName = "FilterCropOutputs";
	string CfolderCreateCommand = "mkdir " + CmyfolderName;
	system(CfolderCreateCommand.c_str());

	filtercropstream << CmyfolderName << "/" << crop << ivar << myzero << segment_no << "_filter" << filter_no << type;
	string Cfilteroutput = filtercropstream.str();
	filtercropstream.str("");
	imwrite(Cfilteroutput, goster_filtre);

	/*clean everything*/
	nindx = 0;
	nIndx.clear();
	for (int ph = 0; ph<pyrHeight; ph++){
		for (int pw = 0; pw<pyrWidth; pw++){
			label[ph][pw] = 0;
			I[ph][pw] = 0;

		}
	}
	for (int ph = 0; ph<pyrHeight; ph++){
		for (int pw = 0; pw<pyrWidth + 1; pw++){

			Q[ph][pw] = 0;

		}
	}

	for (int eqn = 0; eqn<MAX_EXP_NrOf_LABELS; eqn++){ EQ[eqn] = 0; }

	segmentClrs.clear();
	contourArray.clear();
	mychannels.clear();
	fin_img.release();


	myImData.filter.clear();
	myImData.hsv_filter.clear();

	myImData.connComp.clear();
	myImData.original.release();
	myImData.intensity.release();
	myImData.hsvImg.release();
	labelNr = 0;


}


//}
