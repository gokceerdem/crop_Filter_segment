#include "features.h"
#include "SegmentFeatures.h"
#include "SegmentMoments.h"
#include "coor.h"
#include "MyLibConstants.h"
#include "sqrt6.h"

#include <iostream>

#include "ImData.h"

void features(ImData &mid, SegmentFeatures mF, SegmentMoments SegMom, list<SegmentFeatures> featureVector, list<SegmentMoments> momentVector){

	
	//	int segmentNr=0;
	for (int k = 0; k< mid.connComp.size(); k++){

		int cx = 0;
		int cy = 0;

		double filtSum0 = 0;
		double filtSum1 = 0;
		double filtSum2 = 0;
		double filtSum3 = 0;
		double filtSum4 = 0;
		double filtSum5 = 0;
		double filtSum6 = 0;
		double filtSum7 = 0;
		double filtSum8 = 0;
		double filtSum9 = 0;
		double filtSum10 = 0;
		double filtSum11 = 0;
		double filtSum12 = 0;
		double filtSum13 = 0;
		double filtSum14 = 0;
		double filtSum15 = 0;
		double filtSum16 = 0;
		double filtSum17 = 0;
		double filtSum18 = 0;
		double filtSum19 = 0;
		double filtSum20 = 0;
		double filtSum21 = 0;
		double filtSum22 = 0;
		double filtSum23 = 0;
		double filtSum24 = 0;
		double filtSum25 = 0;
		double filtSum26 = 0;
		double filtSum27 = 0;
		double filtSum28 = 0;
		double filtSum29 = 0;
		double filtSum30 = 0;
		double filtSum31 = 0;
		double filtSum32 = 0;
		double filtSum33 = 0;
		double filtSum34 = 0;
		double filtSum35 = 0;
		double filtSum36 = 0;
		double filtSum37 = 0;
		double filtSum38 = 0;
		double filtSum39 = 0;
		double filtSum40 = 0;
		double filtSum41 = 0;
		double filtSum42 = 0;
		double filtSum43 = 0;
		double filtSum44 = 0;
		double filtSum45 = 0;
		double filtSum46 = 0;
		double filtSum47 = 0;
		double filtSum48 = 0;
		double filtSum49 = 0;
		double filtSum50 = 0;
		double filtSum51 = 0;
		double filtSum52 = 0;
		double filtSum53 = 0;
		double filtSum54 = 0; //hue

		int centerxSum = 0;
		int centerySum = 0;

		//second order moments
		double m20 = 0;
		double m02 = 0;
		double m11 = 0;

		double A = mid.connComp.at(k).size(); // Area of the segment = m00
		//	SegMom.m00 = A;

		mF.r = A / (pyrHeight* pyrWidth);

		for (int s = 0; s < mid.connComp.at(k).size(); s++)
		{
			// center of mass calculation --> mass(cx,cy)

			//centerSum = centerSum + memberpixel;
			int memberx = mid.connComp.at(k).at(s).x;
			int membery = mid.connComp.at(k).at(s).y;

			centerxSum = centerxSum + memberx;
			centerySum = centerySum + membery;
			if (s == mid.connComp.at(k).size() - 1)
			{
				//cx and cy calculation
				cx = (centerxSum) / A;
				cy = (centerySum) / A;

			}
		}

		if (cy <= UpObjTh){ mF.segmentLocation = 1; }
		else if (cy>UpObjTh && cy <= MidObjTh){ mF.segmentLocation = 0; }
		else if (cy>MidObjTh){ mF.segmentLocation = -1; }

		for (int s = 0; s < mid.connComp.at(k).size(); s++){

			m20 = m20 + ((mid.connComp.at(k).at(s).x - cx)*(mid.connComp.at(k).at(s).x - cx));
			m02 = m02 + ((mid.connComp.at(k).at(s).y - cy)*(mid.connComp.at(k).at(s).y - cy));
			m11 = m11 + ((mid.connComp.at(k).at(s).x - cx)*(mid.connComp.at(k).at(s).y - cy));


			double dumVal0 = mid.filter.at(0).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum0 = filtSum0 + dumVal0;

			double dumVal1 = mid.filter.at(1).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum1 = filtSum1 + dumVal1;

			double dumVal2 = mid.filter.at(2).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum2 = filtSum2 + dumVal2;

			double dumVal3 = mid.filter.at(3).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum3 = filtSum3 + dumVal3;

			double dumVal4 = mid.filter.at(4).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum4 = filtSum4 + dumVal4;

			double dumVal5 = mid.filter.at(5).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum5 = filtSum5 + dumVal5;

			double	dumVal6 = mid.filter.at(6).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum6 = filtSum6 + dumVal6;

			double	dumVal7 = mid.filter.at(7).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum7 = filtSum7 + dumVal7;

			double	dumVal8 = mid.filter.at(8).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum8 = filtSum8 + dumVal8;

			double	dumVal9 = mid.filter.at(9).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum9 = filtSum9 + dumVal9;

			double	dumVal10 = mid.filter.at(10).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum10 = filtSum10 + dumVal10;

			double	dumVal11 = mid.filter.at(11).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum11 = filtSum11 + dumVal11;

			double	dumVal12 = mid.filter.at(12).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum12 = filtSum12 + dumVal12;

			double	dumVal13 = mid.filter.at(13).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum13 = filtSum13 + dumVal13;

			double	dumVal14 = mid.filter.at(14).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum14 = filtSum14 + dumVal14;

			double	dumVal15 = mid.filter.at(15).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum15 = filtSum15 + dumVal15;

			double	dumVal16 = mid.filter.at(16).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum16 = filtSum16 + dumVal16;

			double	dumVal17 = mid.filter.at(17).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum17 = filtSum17 + dumVal17;

			double	dumVal18 = mid.filter.at(18).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum18 = filtSum18 + dumVal18;

			double	dumVal19 = mid.filter.at(19).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum19 = filtSum19 + dumVal19;

			double	dumVal20 = mid.filter.at(20).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum20 = filtSum20 + dumVal20;

			double	dumVal21 = mid.filter.at(21).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum21 = filtSum21 + dumVal21;

			double	dumVal22 = mid.filter.at(22).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum22 = filtSum22 + dumVal22;

			double	dumVal23 = mid.filter.at(23).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum23 = filtSum23 + dumVal23;

			double	dumVal24 = mid.filter.at(24).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum24 = filtSum24 + dumVal24;

			double	dumVal25 = mid.filter.at(25).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum25 = filtSum25 + dumVal25;

			double	dumVal26 = mid.filter.at(26).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum26 = filtSum26 + dumVal26;

			double	dumVal27 = mid.filter.at(27).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum27 = filtSum27 + dumVal27;

			double	dumVal28 = mid.filter.at(28).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum28 = filtSum28 + dumVal28;

			double	dumVal29 = mid.filter.at(29).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum29 = filtSum29 + dumVal29;

			double	dumVal30 = mid.filter.at(30).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum30 = filtSum30 + dumVal30;

			double	dumVal31 = mid.filter.at(31).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum31 = filtSum31 + dumVal31;

			double	dumVal32 = mid.filter.at(32).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum32 = filtSum32 + dumVal32;

			double	dumVal33 = mid.filter.at(33).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum33 = filtSum33 + dumVal33;

			double	dumVal34 = mid.filter.at(34).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum34 = filtSum34 + dumVal34;

			double	dumVal35 = mid.filter.at(35).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum35 = filtSum35 + dumVal35;

			double	dumVal36 = mid.filter.at(36).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum36 = filtSum36 + dumVal36;

			double	dumVal37 = mid.filter.at(37).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum37 = filtSum37 + dumVal37;

			double	dumVal38 = mid.filter.at(38).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum38 = filtSum38 + dumVal38;

			double	dumVal39 = mid.filter.at(39).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum39 = filtSum39 + dumVal39;

			double	dumVal40 = mid.filter.at(40).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum40 = filtSum40 + dumVal40;

			double	dumVal41 = mid.filter.at(41).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum41 = filtSum41 + dumVal41;

			double	dumVal42 = mid.filter.at(42).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum42 = filtSum42 + dumVal42;

			double	dumVal43 = mid.filter.at(43).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum43 = filtSum43 + dumVal43;

			double	dumVal44 = mid.filter.at(44).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum44 = filtSum44 + dumVal44;

			double	dumVal45 = mid.filter.at(45).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum45 = filtSum45 + dumVal45;

			double	dumVal46 = mid.filter.at(46).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum46 = filtSum46 + dumVal46;

			double	dumVal47 = mid.filter.at(47).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum47 = filtSum47 + dumVal47;

			double	dumVal48 = mid.filter.at(48).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum48 = filtSum48 + dumVal48;

			double	dumVal49 = mid.filter.at(49).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum49 = filtSum49 + dumVal49;

			double	dumVal50 = mid.filter.at(50).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum50 = filtSum50 + dumVal50;

			double	dumVal51 = mid.filter.at(51).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum51 = filtSum51 + dumVal51;

			double	dumVal52 = mid.filter.at(52).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum52 = filtSum52 + dumVal52;

			double	dumVal53 = mid.filter.at(53).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum53 = filtSum53 + dumVal53;

			double	dumVal54 = mid.filter.at(54).at<uchar>(mid.connComp.at(k).at(s).y, mid.connComp.at(k).at(s).x);
			filtSum54 = filtSum54 + dumVal54;

			if (s == mid.connComp.at(k).size() - 1)
			{
				// normalized feature results
				double f0 = filtSum0 / A;
				f0 = f0 / 255;
				mF.AverageFilterOutputs.push_back(f0);

				double f1 = filtSum1 / A;
				f1 = f1 / 255;
				mF.AverageFilterOutputs.push_back(f1);

				double f2 = filtSum2 / A;
				f2 = f2 / 255;
				mF.AverageFilterOutputs.push_back(f2);

				double f3 = filtSum3 / A;
				f3 = f3 / 255;
				mF.AverageFilterOutputs.push_back(f3);

				double f4 = filtSum4 / A;
				f4 = f4 / 255;
				mF.AverageFilterOutputs.push_back(f4);

				double f5 = filtSum5 / A;
				f5 = f5 / 255;
				mF.AverageFilterOutputs.push_back(f5);

				double	f6 = filtSum6 / A;
				f6 = f6 / 255;
				mF.AverageFilterOutputs.push_back(f6);

				double	f7 = filtSum7 / A;
				f7 = f7 / 255;
				mF.AverageFilterOutputs.push_back(f7);

				double	f8 = filtSum8 / A;
				f8 = f8 / 255;
				mF.AverageFilterOutputs.push_back(f8);

				double	f9 = filtSum9 / A;
				f9 = f9 / 255;
				mF.AverageFilterOutputs.push_back(f9);

				double	f10 = filtSum10 / A;
				f10 = f10 / 255;
				mF.AverageFilterOutputs.push_back(f10);

				double	f11 = filtSum11 / A;
				f11 = f11 / 255;
				mF.AverageFilterOutputs.push_back(f11);

				double	f12 = filtSum12 / A;
				f12 = f12 / 255;
				mF.AverageFilterOutputs.push_back(f12);

				double	f13 = filtSum13 / A;
				f13 = f13 / 255;
				mF.AverageFilterOutputs.push_back(f13);

				double	f14 = filtSum14 / A;
				f14 = f14 / 255;
				mF.AverageFilterOutputs.push_back(f14);

				double	f15 = filtSum15 / A;
				f15 = f15 / 255;
				mF.AverageFilterOutputs.push_back(f15);

				double	f16 = filtSum16 / A;
				f16 = f16 / 255;
				mF.AverageFilterOutputs.push_back(f16);

				double	f17 = filtSum17 / A;
				f17 = f17 / 255;
				mF.AverageFilterOutputs.push_back(f17);

				double	f18 = filtSum18 / A;
				f18 = f18 / 255;
				mF.AverageFilterOutputs.push_back(f18);

				double	f19 = filtSum19 / A;
				f19 = f19 / 255;
				mF.AverageFilterOutputs.push_back(f19);

				double	f20 = filtSum20 / A;
				f20 = f20 / 255;
				mF.AverageFilterOutputs.push_back(f20);

				double	f21 = filtSum21 / A;
				f21 = f21 / 255;
				mF.AverageFilterOutputs.push_back(f21);

				double	f22 = filtSum22 / A;
				f22 = f22 / 255;
				mF.AverageFilterOutputs.push_back(f22);

				double	f23 = filtSum23 / A;
				f23 = f23 / 255;
				mF.AverageFilterOutputs.push_back(f23);

				double	f24 = filtSum24 / A;
				f24 = f24 / 255;
				mF.AverageFilterOutputs.push_back(f24);

				double	f25 = filtSum25 / A;
				f25 = f25 / 255;
				mF.AverageFilterOutputs.push_back(f25);

				double	f26 = filtSum26 / A;
				f26 = f26 / 255;
				mF.AverageFilterOutputs.push_back(f26);

				double	f27 = filtSum27 / A;
				f27 = f27 / 255;
				mF.AverageFilterOutputs.push_back(f27);

				double	f28 = filtSum28 / A;
				f28 = f28 / 255;
				mF.AverageFilterOutputs.push_back(f28);

				double	f29 = filtSum29 / A;
				f29 = f29 / 255;
				mF.AverageFilterOutputs.push_back(f29);

				double	f30 = filtSum30 / A;
				f30 = f30 / 255;
				mF.AverageFilterOutputs.push_back(f30);

				double	f31 = filtSum31 / A;
				f31 = f31 / 255;
				mF.AverageFilterOutputs.push_back(f31);

				double	f32 = filtSum32 / A;
				f32 = f32 / 255;
				mF.AverageFilterOutputs.push_back(f32);

				double	f33 = filtSum33 / A;
				f33 = f33 / 255;
				mF.AverageFilterOutputs.push_back(f33);

				double	f34 = filtSum34 / A;
				f34 = f34 / 255;
				mF.AverageFilterOutputs.push_back(f34);

				double	f35 = filtSum35 / A;
				f35 = f35 / 255;
				mF.AverageFilterOutputs.push_back(f35);

				double	f36 = filtSum36 / A;
				f36 = f36 / 255;
				mF.AverageFilterOutputs.push_back(f36);

				double	f37 = filtSum37 / A;
				f37 = f37 / 255;
				mF.AverageFilterOutputs.push_back(f37);

				double	f38 = filtSum38 / A;
				f38 = f38 / 255;
				mF.AverageFilterOutputs.push_back(f38);

				double	f39 = filtSum39 / A;
				f39 = f39 / 255;
				mF.AverageFilterOutputs.push_back(f39);

				double	f40 = filtSum40 / A;
				f40 = f40 / 255;
				mF.AverageFilterOutputs.push_back(f40);

				double	f41 = filtSum41 / A;
				f41 = f41 / 255;
				mF.AverageFilterOutputs.push_back(f41);

				double	f42 = filtSum42 / A;
				f42 = f42 / 255;
				mF.AverageFilterOutputs.push_back(f42);

				double	f43 = filtSum43 / A;
				f43 = f43 / 255;
				mF.AverageFilterOutputs.push_back(f43);

				double	f44 = filtSum44 / A;
				f44 = f44 / 255;
				mF.AverageFilterOutputs.push_back(f44);

				double	f45 = filtSum45 / A;
				f45 = f45 / 255;
				mF.AverageFilterOutputs.push_back(f45);

				double	f46 = filtSum46 / A;
				f46 = f46 / 255;
				mF.AverageFilterOutputs.push_back(f46);

				double	f47 = filtSum47 / A;
				f47 = f47 / 255;
				mF.AverageFilterOutputs.push_back(f47);

				double	f48 = filtSum48 / A;
				f48 = f48 / 255;
				mF.AverageFilterOutputs.push_back(f48);

				double	f49 = filtSum49 / A;
				f49 = f49 / 255;
				mF.AverageFilterOutputs.push_back(f49);

				double	f50 = filtSum50 / A;
				f50 = f50 / 255;
				mF.AverageFilterOutputs.push_back(f50);

				double	f51 = filtSum51 / A;
				f51 = f51 / 255;
				mF.AverageFilterOutputs.push_back(f51);

				double	f52 = filtSum52 / A;
				f52 = f52 / 255;
				mF.AverageFilterOutputs.push_back(f52);

				double	f53 = filtSum53 / A;
				f53 = f53 / 255;
				mF.AverageFilterOutputs.push_back(f53);


				//hue normalization is different
				double f54 = filtSum54 / A;
				f54 = f54 / 180;
				mF.AverageFilterOutputs.push_back(f54);

				//mF.mass = Point(cx,cy);

				//calculation of major and minor axis

				m20 = m20 / A; SegMom.m20 = m20;
				m02 = m02 / A; SegMom.m02 = m02;
				m11 = m11 / A; SegMom.m11 = m11;


				double k4 = sqrt6(abs((4 * (m11)*(m11)) - ((m20 - m02)*(m20 - m02))));

				double lambda1 = sqrt6(((m20 + m02) / 2) + k4); //major
				lambda1 = lambda1 / hypo;	// = q2 (normalized major length)
				double lambda2 = sqrt6(abs(((m20 + m02) / 2) - k4)); //minor
				lambda2 = lambda2 / hypo;	// =  q3 (normalized minor length)


				mF.elongation.push_back(lambda1);
				mF.elongation.push_back(lambda2);
				//eccentricity = abs((((m20 - m02)*(m20 - m02))-(4*(m11)*(m11)))/((m20 + m02)*(m20 + m02)));

				featureVector.push_back(mF);
				//	momentVector.push_back(SegMom);

				mF.AverageFilterOutputs.clear();
				mF.elongation.clear();

				//				segmentNr = segmentNr+1;
			}
		}
	}
}