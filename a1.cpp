#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <DrawText.h>

//Cmmented for check
using namespace std;
////Lets see
// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) s. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them 
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);
      
    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}

// DetectedSymbol class may be helpful!
//  Feel free to modify.
//
typedef enum {NOTEHEAD=0, QUARTERREST=1, EIGHTHREST=2} Type;
class DetectedSymbol {
public:
  int row, col, width, height;
  Type type;
  char pitch;
  double confidence;
};

// Function that outputs the ascii detection output file
void  write_detection_txt(const string &filename, const vector<struct DetectedSymbol> &symbols)
{
  ofstream ofs(filename.c_str());

  for(int i=0; i<symbols.size(); i++)
    {
      const DetectedSymbol &s = symbols[i];
      ofs << s.row << " " << s.col << " " << s.width << " " << s.height << " ";
      if(s.type == NOTEHEAD)
	ofs << "filled_note " << s.pitch;
      else if(s.type == EIGHTHREST)
	ofs << "eighth_rest _";
      else 
	ofs << "quarter_rest _";
      ofs << " " << s.confidence << endl;
    }
}

// Function that outputs a visualization of detected symbols
void  write_detection_image(const string &filename, const vector<DetectedSymbol> &symbols, const SDoublePlane &input)
{
  SDoublePlane output_planes[3];
  for(int i=0; i<3; i++)
    output_planes[i] = input;

  for(int i=0; i<symbols.size(); i++)
    {
      const DetectedSymbol &s = symbols[i];

      overlay_rectangle(output_planes[s.type], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 255, 2);
      overlay_rectangle(output_planes[(s.type+1) % 3], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 0, 2);
      overlay_rectangle(output_planes[(s.type+2) % 3], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 0, 2);

      if(s.type == NOTEHEAD)
	{
	  char str[] = {s.pitch, 0};
	  draw_text(output_planes[0], str, s.row, s.col+s.width+1, 0, 2);
	  draw_text(output_planes[1], str, s.row, s.col+s.width+1, 0, 2);
	  draw_text(output_planes[2], str, s.row, s.col+s.width+1, 0, 2);
	}
    }

  SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}



// The rest of these functions are incomplete. These are just suggestions to 
// get you started -- feel free to add extra functions, change function
// parameters, etc.

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
    SDoublePlane output(input.rows(), input.cols());
    SDoublePlane sepin(input.rows()+(2*(col_filter.rows()/2)), input.cols()+(2*(row_filter.cols()/2)));
    // Convolution code here
    
    for(int x=0;x<input.rows();x++)
    {
        for(int y=0;y<input.cols();y++)          {
            sepin[x+(col_filter.rows()/2)][y+(row_filter.cols()/2)]=input[x][y];
        }
    }
    
    for(int i=0;i<output.rows();i++){
        for(int j=0;j<output.cols();j++){
            for(int k=0;k<col_filter.rows();k++){
                output[i][j]+=sepin[i+k][j+1]*col_filter[col_filter.rows()-1-k][0];
            }
        }
    }
    
    
    for(int x=0;x<sepin.rows();x++)
    {
        for(int y=0;y<sepin.cols();y++)          {
            sepin[x][y]=0;
        }
    }
    //SDoublePlane sepin2(input.rows(), input.cols()+(2*(row_filter.cols()/2)));
    for(int x=0;x<input.rows();x++)
    {
        for(int y=0;y<input.cols();y++)          {
            sepin[x+(col_filter.rows()/2)][y+(row_filter.cols()/2)]=output[x][y];
            output[x][y]=0;
        }
    }
    
    for(int i=0;i<output.rows();i++){
        for(int j=0;j<output.cols();j++){
            for(int k=0;k<row_filter.cols();k++){
                output[i][j]+=sepin[i+1][j+k]*row_filter[0][row_filter.cols()-1-k];
            }
        }
    }
    
    return output;
}

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
    int fci = filter.rows()/2;
    
    int fcj = filter.cols()/2;
    
    SDoublePlane in(input.rows()+(2*fci), input.cols()+(2*fcj));
    
    for(int x=0;x<input.rows();x++)
    {
        for(int y=0;y<input.cols();y++)          {
            in[x+fci][y+fcj]=input[x][y];
        }
    }
    
    SDoublePlane output(input.rows(), input.cols());
    
    for(int x=0;x<output.rows();x++)
    {
        for(int y=0;y<output.cols();y++)          {
            
            for(int i=0;i<filter.rows();i++)          {
                for(int j=0;j<filter.cols();j++)          {
                    output[x][y]+=in[x+i][y+j]*filter[filter.rows()-1-i][filter.cols()-1-j];
                }
                
            }
            //cout<<x<<" "<<y<<" "<<endl;

            
        }
        //cout<<endl;
    }

//delete in;
  //  free (in);
  // Convolution code here
    
    
  
  return output;
}


// Apply a sobel operator to an image, returns the result
// 
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
 // SDoublePlane output(input.rows(), input.cols());

  // Implement a sobel gradient estimation filter with 1-d filters
    
    SDoublePlane sobely(3,3);
    //sobel={{-1,0,1},{-2,0,2},{-1,0,1}};
    sobely[0][0]=-1;
    sobely[0][2]=1;
    sobely[1][0]=-2;
    sobely[1][2]=2;
    sobely[2][0]=-1;
    sobely[2][2]=1;
    SDoublePlane outputy=convolve_general(input, sobely);
    SDoublePlane sobelx(3,3);
    //sobel={{-1,-2,-1},{0,0,0},{-1,-2,-1}};
    sobelx[0][0]=-1;
    sobelx[0][1]=-2;
    sobelx[0][2]=-1;
    sobelx[2][0]=1;
    sobelx[2][1]=2;
    sobelx[2][2]=1;
    
    SDoublePlane outputx=convolve_general(input, sobelx);
    SDoublePlane output(input.rows(),input.cols());
    double max=0;
    for(int i=0;i<input.rows();i++)
        for(int j=0;j<input.cols();j++)
        {output[i][j]=sqrt((outputx[i][j]*outputx[i][j])+(outputy[i][j]*outputy[i][j]));
            if(max<output[i][j]){
                max=output[i][j];
            }
        }
    for(int i=0;i<input.rows();i++)
        for(int j=0;j<input.cols();j++)
        {output[i][j]=output[i][j]*255/max;
            
        }
  return output;
}

// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(const SDoublePlane &input, double thresh=0)
{
  SDoublePlane output(input.rows(), input.cols());

  // Implement an edge detector of your choice, e.g.
  // use your sobel gradient operator to compute the gradient magnitude and threshold
  
  return output;
}

SDoublePlane to_binary_image(const SDoublePlane &input)
{
    SDoublePlane output(input.rows(), input.cols());
    for(int x=0;x<input.rows();x++)
        for(int y=0;y<input.cols();y++)
        {
            if(input[x][y]<128)
            output[x][y]=0;
            else
                output[x][y]=1;

            
        }
    return output;
}

SDoublePlane to_binary_image_minus1(const SDoublePlane &input)
{
    SDoublePlane output(input.rows(), input.cols());
    for(int x=0;x<input.rows();x++)
        for(int y=0;y<input.cols();y++)
        {
            
            if(input[x][y]<128)
                output[x][y]=1;
            else
                output[x][y]=0;
        }
    return output;
}

SDoublePlane matching(const SDoublePlane &input_image,const SDoublePlane &template1_image)
{
    int fci=template1_image.rows()/2;
    int fcj=template1_image.cols()/2;
    
    SDoublePlane match(input_image.rows(),input_image.cols());
    
    double totalpixelsintempnorm255=255.0/(template1_image.rows()*template1_image.cols());
    
    for(int x=fci;x<input_image.rows()-fci;x++)
        for(int y=fcj;y<input_image.cols()-fcj;y++)
        {
            for(int i=0;i<template1_image.rows();i++)
                    for(int j=0;j<template1_image.cols();j++)
                    {if(input_image[x-fci+i][y-fcj+j]<128 && template1_image[i][j]<128)
                        match[x][y]+=totalpixelsintempnorm255;
                     else   if(input_image[x-fci+i][y-fcj+j]>=128 && template1_image[i][j]>=128)
                            match[x][y]+=totalpixelsintempnorm255;}
        }
    return match;
 }
//
// This main file just outputs a few test images. You'll want to change it to do 
//  something more interesting!
//

SDoublePlane matching_alt(const SDoublePlane &input_image,const SDoublePlane &template1_image){
    SDoublePlane match(input_image.rows(),input_image.cols());
    //SDoublePlane gamma_mat(input_image.rows(),input_image.cols());
    SDoublePlane dist_mat(input_image.rows(),input_image.cols());
    SDoublePlane temp_img(template1_image.rows(),template1_image.cols());
    //gamma_mat = input_image;
    //gamma_mat = gamma_func(gamma_mat);
    //SDoublePlane input_binary(input_image.rows(),input_image.cols());
    //SDoublePlane temp_binary(template1_image.rows(),template1_image.cols());
    SDoublePlane input_binary = to_binary_image(input_image);
    SDoublePlane temp_binary = to_binary_image(template1_image);
    //Initialization done
    for(int i=0;i<input_binary.rows();i++){
        for(int j=0;j<input_binary.cols();j++){
            //match[i][j]=1;
            if(input_binary[i][j]==0){
                dist_mat[i][j]=1;
            }
            else{
                dist_mat[i][j]=0;
            }
        }
    }
    //dist_mat = gamma_mat;
    
    for(int i=0;i<template1_image.rows();i++){
        for(int j=0;j<template1_image.cols();j++){
            if(template1_image[i][j]!=0){
                temp_img[i][j]=1;
            }
            else{
                temp_img[i][j]=0;
            }
        }
    }
    
    //quadratic complexity to find the smallest distance from the edge(kernel from Cornell)
    for(int i=1;i<dist_mat.rows();i++){
        for(int j=1;j<dist_mat.cols();j++){
            //dist_mat[i][j]=min(min(min(dist_mat[i][j],dist_mat[i-1][j]+1),dist_mat[i][j-1]+1),dist_mat[i-1][j-1]+sqrt(2));
            dist_mat[i][j]=min(min(dist_mat[i][j],dist_mat[i-1][j]+1),dist_mat[i][j-1]+1);
        }
    }
    
    
    for(int i=dist_mat.rows()-2;i>=0;i--){
        for(int j=dist_mat.cols()-2;j>=0;j--){
            //dist_mat[i][j]=min(min(min(dist_mat[i][j],dist_mat[i+1][j]+1),dist_mat[i][j+1]+1),dist_mat[i+1][j+1]+sqrt(2));
            dist_mat[i][j]=min(min(dist_mat[i][j],dist_mat[i+1][j]+1),dist_mat[i][j+1]+1);
        }
    }
    dist_mat[0][dist_mat.cols()-1] = min(min(dist_mat[0][dist_mat.cols()-1],dist_mat[0][dist_mat.cols()-2]+1),dist_mat[1][dist_mat.cols()-1]+1);
    dist_mat[dist_mat.rows()-1][0]= min(min(dist_mat[dist_mat.rows()-1][0],dist_mat[dist_mat.rows()-2][0]+1),dist_mat[dist_mat.rows()-1][1]+1);
    
    int fci = temp_binary.rows()/2;
    int fcj = temp_binary.cols()/2;
    double outmax = 0;
    for(int i=fci;i<match.rows()-fci;i++){
        for(int j=fcj;j<match.cols()-fcj;j++){
            
            for(int k=0;k<temp_binary.rows();k++){
                for(int l=0;l<temp_binary.cols();l++){
                    
                    match[i][j]+=(temp_binary[k][l]*dist_mat[i-(fci)+k][j-(fcj)+l]);
                    
                }
            }
            

            if(outmax<match[i][j]){
                outmax = match[i][j];
            }
        }
    }
    //cout<<outmax;
    //cout<<outmax;
    SDoublePlane matchfin((match.rows()-(2*fci)),(match.cols()-(2*fcj)));
    //cout<<"hey";
    double max=0;
    for(int i=fci;i<match.rows()-fci;i++){
        for(int j=fcj;j<match.cols()-fcj;j++){
            
            //cout<<i<<"hey"<<j<<endl;
            
            
            matchfin[i-fci][j-fci]=match[i][j];
            if(max<match[i][j]){
                max = match[i][j];
            }
        }
    }
    //cout<<max<<endl;
    //cout<<outmax<<endl;
    	for(int i=0;i<matchfin.rows();i++){
            for(int j=0;j<matchfin.cols();j++){
                    matchfin[i][j]=255-(matchfin[i][j]*255/max);}
     }
    
    return matchfin;
}
/*SDoublePlane HT(const SDoublePlane &input_image) {
    SDoublePlane ip = to_binary_image(input_image);
    
    vector<DetectedSymbol> workaround1;
    write_detection_image("ht1.png",workaround1,ip);
    int space=ip.rows()/10;
    cout<<"Space "<<space<<endl;
    int r=ip.rows();
    cout<<"r "<<r<<endl;

    SDoublePlane HA(r-space,space);
    cout<<"HArows "<<HA.rows()<<" HAcols "<<HA.cols()<<endl;
    double max=0;
    int h = HA.rows();
    int w = ip.cols();
    cout<<ip.rows()<<" "<<h<<endl;
    cout<<ip.cols()<<" "<<w<<endl;
    for (int v = 0; v < HA.rows(); v++) {
        cout<<"hello"<<endl;
        for (int u = 0; u < ip.cols(); u++) {
            if (ip[u][v] > 0) {				//edge pixels should be given a value of 1, to distinguish.
                //	doThis(u, v);
                //}
               cout<<"hi"<<endl;
                //int x = u - xCtr, y = v - yCtr;
                for (int i = 1; i < HA.cols(); i++) {
                    cout<<"hey"<<i<<" "<<u<<" "<<v<<endl;
                    cout<<ip.rows()<<endl;
                    cout<<HA.rows()<<" "<<v<<endl<<HA.cols()<<" "<<i<<endl;
                    //cout<<i<<endl;
                    //double theta = dAng * ia;
                    
                    if(v+4*i < ip.rows()) {
                        cout<<"yo";
                    double j = ip[u][v] + ip[u][v+i] + ip[u][v+2*i] + ip[u][v+3*i] + ip[u][v+4*i];//cRad + (int) round((x*cos(theta) + y*sin(theta)) / dRad);
                        cout<<"yo";

                    if (j > 0) {
                    
                        HA[v][i]++;
                        if (max<HA[v][i])
                            max=HA[v][i];
                    }
                    }
                }//cout<< HA[ia][ir];
                
            }
            //cout<<u<<endl;
        }	
    }
    cout<<max;
    for (int v = 0; v < HA.rows(); v++) {
        for (int u = 0; u < HA.cols(); u++) {
            HA[v][u]=HA[v][u]*255/max;
        }
    }
    return HA;
}
*/

SDoublePlane HT(const SDoublePlane &input_image) {
    SDoublePlane ip = to_binary_image(input_image);
    
    vector<DetectedSymbol> workaround1;
    write_detection_image("ht1.png",workaround1,ip);
    //int space=ip.rows()/10;
    //cout<<"Space "<<space<<endl;
    int r=ip.rows();
    
    
    SDoublePlane HA(r,1);
    int start[r];
    for(int i=0;i<r;i++)
    {
        start[i]=0;
    }
    double max=0;
    for(int i=1;i<ip.rows()-1;i++)
    {
        for(int j=0;j<ip.cols();j++)
        {if(ip[i][j]>0)
        {
            start[i]=start[i]+1;
            HA[i][0]+=1;
            if(HA[i][0]>max)
                max=HA[i][0];
        }
        }

    }
    for(int i=0;i<r;i++)
    {
        HA[i][0]=HA[i][0]*255/max;
    }
      return HA;
}

int main(int argc, char *argv[])
{
  if(!(argc == 2))
    {
      cerr << "usage: " << argv[0] << " input_image" << endl;
      return 1;
    }
    

  string input_filename(argv[1]);
  SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());
    SDoublePlane template1_image= SImageIO::read_png_file("t1.png");
    SDoublePlane template2_image= SImageIO::read_png_file("template2.png");
    SDoublePlane template3_image= SImageIO::read_png_file("template3.png");
    
    
    
   
  vector<DetectedSymbol> symbols;
    vector<DetectedSymbol> symbols5;
    vector<DetectedSymbol> workaround;
    SDoublePlane inputedge=sobel_gradient_filter(input_image,true);

    
  //  cout<<"HT";
    SDoublePlane h=HT(inputedge);
    //write_detection_image("ht.png", workaround, h);
    SDoublePlane HTop(input_image.rows(),input_image.cols());
    for(int i=1;i<h.rows()-1;i++)
    {
        //cout<<h[i][0]<<endl;
        for(int j=0;j<HTop.cols();j++)
            if(h[i][0]>200)
            {
                bool falsepos=false;

                for(int k=1;k<5;k++)
                    if(h[i+k][0]>200)
                    {falsepos=true;
                        break;
                    }
                if(!falsepos)
                {HTop[i][j]=255;
                    h[i][0]=255;
                }
                else
                    h[i][0]=0;
            }
        else
            h[i][0]=0;
    }
    
    write_detection_image("staves.png", workaround, HTop);
   // cout<<"HTend"<<endl;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
 
    
    SDoublePlane tempedge1=sobel_gradient_filter(template1_image,true);
    SDoublePlane tempedge2=sobel_gradient_filter(template2_image,true);

    SDoublePlane tempedge3=sobel_gradient_filter(template3_image,true);


    write_detection_image("edges.png",symbols,inputedge);
   // write_detection_image("tempedge.png",symbols,tempedge);
    
    //SDoublePlane H=HT(inputedge);
    //write_detection_image("ht.png",workaround,H);
    
    
 
    

    
    
     SDoublePlane matchdraw=matching(input_image,template1_image);
     write_detection_image("scores41.png",workaround,matchdraw);

    
    
    
    int dis=1000;
    int first[4]={0,0,0,0};
    bool found=false;
    int k=0;
    for(int i=0;i<h.rows();i++)
    {
        //cout<<first[k]<<" "<<dis<<" "<<i<<endl;
        if(h[i][0]==255)
        {
            if(!found){
                first[k]=i;
                //cout<<"hey"<<first[k];
                k++;
                
                found=true;
                if(k>3)
                    break;
                
            }
            else
            {
                int temp=i-first[k-1];
                //cout<<"temp"<<temp<<endl;
                if(temp<dis)
                    dis=temp;
                found=false;
                if(first[k-1]+(dis*5)<h.rows())
                i=first[k-1]+(dis*5);
                
            }
        }
    }
    //cout<<first[0]<<first[1]<<first[2]<<first[3]<<endl;
    ////cout<<dis<<endl;
    
    for(int i=0;i<matchdraw.rows();i++)    {
        for(int j=0;j<matchdraw.cols();j++)    {
            bool falsepos=false;

            if(matchdraw[i][j]>200)
            {
                for(int x=i-template1_image.rows()/2;x<i+template1_image.rows()/2 && !falsepos;x++)
                    for(int y=j-template1_image.cols();y<j+template1_image.cols()/2;y++)
                        if(matchdraw[x][y]>matchdraw[i][j])
                        {falsepos=true;
                           // cout<<"heyin"<<falsepos<<endl;
                            break;
                        }
             
                //cout<<"heyout"<<falsepos<<endl;

      if(!falsepos)
      {

          DetectedSymbol s;
          int ffirst=0;
          int t=1000;
          for(int a=0;a<4;a++)
          {
              if(abs(first[a]-i)<t)
              {
                  t=abs(first[a]-i);
                  ffirst=first[a];
              }
          }
         // cout<<ffirst;
          int g=1000;
          int pos;
          for(int i=0;i<10;i++)
          {
              if(abs(ffirst+(i*dis/2))<t)
              {g=abs(ffirst+(i*dis/2));
                  pos=i;}
              
              
          }
          int note;
          if(pos==0)
              note=5;
          if(pos==1)
              note=4;
          if(pos==2)
              note=3;
          if(pos==3)
              note=2;
          if(pos==4)
              note=1;
          if(pos==5)
              note=0;
          if(pos==6)
              note=6;
          if(pos==7)
              note=5;
          if(pos==8)
              note=4;
          if(pos==9)
              note=3;
          
      s.row = i-template1_image.rows()/2;
      s.col = j-template1_image.cols()/2;
      s.width = template1_image.cols();
      s.height = template1_image.rows();
      s.type = (Type) 0;
      s.confidence = matchdraw[i][j];
      s.pitch = note + 'A';//"Change this as per note"
          symbols.push_back(s);}}
        }
    }


    
    //delete(template1_image);
    //delete(matchdraw1);
  
     matchdraw=matching(input_image,template2_image);
    write_detection_image("scores42.png",workaround,matchdraw);

    for(int i=0;i<matchdraw.rows();i++)    {
        for(int j=0;j<matchdraw.cols();j++)    {
            bool falsepos=false;
            
            if(matchdraw[i][j]>220)
            {
                for(int x=i-template2_image.rows()/2;x<i+template2_image.rows()/2 && !falsepos;x++)
                    for(int y=j-template2_image.cols();y<j+template2_image.cols()/2;y++)
                        if(matchdraw[x][y]>matchdraw[i][j])
                        {falsepos=true;
                            // cout<<"heyin"<<falsepos<<endl;
                            break;
                        }
                
                //cout<<"heyout"<<falsepos<<endl;
                
                if(!falsepos)
                {
                    
                    DetectedSymbol s;
                    s.row = i-template2_image.rows()/2-1;
                    s.col = j-template2_image.cols()/2-1;
                    s.width = template2_image.cols()+2;
                    s.height = template2_image.rows()+2;
                    s.type = (Type) 1;
                    s.confidence = matchdraw[i][j];
                    s.pitch = (rand() % 7) + 'A';//"Change this as per note"
                    symbols.push_back(s);}}
        }
    }

    
    
     matchdraw=matching(input_image,template3_image);
    write_detection_image("scores43.png",workaround,matchdraw);
//
    for(int i=0;i<matchdraw.rows();i++)    {
        for(int j=0;j<matchdraw.cols();j++)    {
            bool falsepos=false;
            
            if(matchdraw[i][j]>220)
            {
                for(int x=i-template3_image.rows()/2;x<i+template3_image.rows()/2 && !falsepos;x++)
                    for(int y=j-template3_image.cols();y<j+template3_image.cols()/2;y++)
                        if(matchdraw[x][y]>matchdraw[i][j])
                        {falsepos=true;
                            break;
                        }
                
                
                if(!falsepos)
                {
                    
                    DetectedSymbol s;
                    s.row = i-template3_image.rows()/2;
                    s.col = j-template3_image.cols()/2;
                    s.width = template3_image.cols();
                    s.height = template3_image.rows();
                    s.type = (Type) 2;
                    s.confidence = matchdraw[i][j];
                    s.pitch = (rand() % 7) + 'A';//"Change this as per note"
                    symbols.push_back(s);}}
        }
    }
 //   delete(template3_image);
   // delete(matchdraw3);

  write_detection_txt("detected7.txt", symbols);
    write_detection_txt("detected4.txt", symbols);

    write_detection_image("detected4.png", symbols, input_image);
    write_detection_image("detected7.png", symbols, input_image);
    
    SDoublePlane matchdraw1=matching_alt(inputedge,tempedge1);
  // write_detection_image("scores51.png",workaround,matchdraw1);
    
    for(int i=0;i<matchdraw1.rows();i++)    {
        for(int j=0;j<matchdraw1.cols();j++)    {
            bool falsepos=false;
            
            if(matchdraw1[i][j]>150)
            {
                for(int x=i-tempedge1.rows()/2;x<i+tempedge1.rows()/2 && !falsepos && x>0 && x<matchdraw1.rows();x++)
                    for(int y=j-tempedge1.cols();y<j+tempedge1.cols()/2 && y>0 && y<matchdraw1.cols();y++)
                        if(matchdraw1[x][y]>matchdraw1[i][j])
                        {falsepos=true;
                            // cout<<"heyin"<<falsepos<<endl;
                            break;
                        }
                
                //cout<<"heyout"<<falsepos<<endl;
                
                if(!falsepos)
                {
                    
                    DetectedSymbol s;
                    s.row = i;
                    s.col = j;
                    s.width = tempedge1.cols();
                    s.height = tempedge1.rows();
                    s.type = (Type) 0;
                    s.confidence = matchdraw1[i][j];
                    s.pitch = (rand() % 7) + 'A';//"Change this as per note"
                    symbols5.push_back(s);}}
        }
    }
     matchdraw1=matching_alt(inputedge,tempedge2);
 //   write_detection_image("scores52.png",workaround,matchdraw1);
    
    for(int i=0;i<matchdraw1.rows();i++)    {
        for(int j=0;j<matchdraw1.cols();j++)    {
            bool falsepos=false;
            
            if(matchdraw1[i][j]>150)
            {
                for(int x=i-tempedge2.rows()/2;x<i+tempedge2.rows()/2 && !falsepos && x>0 && x<matchdraw1.rows();x++)
                    for(int y=j-tempedge2.cols();y<j+tempedge2.cols()/2 && y>0 && y<matchdraw1.cols();y++)
                        if(matchdraw1[x][y]>matchdraw1[i][j])
                        {falsepos=true;
                            break;
                        }
                
                
                if(!falsepos)
                {
                    
                    DetectedSymbol s;
                    s.row = i;
                    s.col = j;
                    s.width = tempedge2.cols();
                    s.height = tempedge2.rows();
                    s.type = (Type) 1;
                    s.confidence = matchdraw1[i][j];
                    s.pitch = (rand() % 7) + 'A';
                    symbols5.push_back(s);}}
        }
    }
     matchdraw1=matching_alt(inputedge,tempedge3);
  // write_detection_image("scores53.png",workaround,matchdraw1);
    
    for(int i=0;i<matchdraw1.rows();i++)    {
        for(int j=0;j<matchdraw1.cols();j++)    {
            bool falsepos=false;
            
            if(matchdraw1[i][j]>150)
            {
                for(int x=i-tempedge3.rows()/2;x<i+tempedge3.rows()/2 && !falsepos && x>0 && x<matchdraw1.rows();x++)
                    for(int y=j-tempedge3.cols();y<j+tempedge3.cols()/2 && y>0 && y<matchdraw1.cols();y++)
                        if(matchdraw1[x][y]>matchdraw1[i][j])
                        {falsepos=true;
                            // cout<<"heyin"<<falsepos<<endl;
                            break;
                        }
                
                //cout<<"heyout"<<falsepos<<endl;
                
                if(!falsepos)
                {
                    
                    DetectedSymbol s;
                    s.row = i;
                    s.col = j;
                    s.width = tempedge3.cols();
                    s.height = tempedge3.rows();
                    s.type = (Type) 2;
                    s.confidence = matchdraw1[i][j];
                    s.pitch = (rand() % 7) + 'A';//"Change this as per note"
                    symbols5.push_back(s);}}
        }
    }
    write_detection_txt("detected5.txt", symbols5);
    
    write_detection_image("detected5.png", symbols5, input_image);
    

}
