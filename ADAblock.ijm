// ImageJ BCP Structure Analysis Framework  

//{{{ collapse all...

// ABOUT THE CODE

	/* ADAblock: Automated Defect Analysis of Block Copolymers
	 *  
	 *   ImageJ Macro Code Implementation of Defect Analysis Algorithm
	 *   Version: 0.50i  Date: 2015.01.30  Author: Jeffrey N. Murphy
	 *   Updated versions available at: https://github.com/MurphysLab/ADAblock
	 *   
	 * The algorithm is described in the following paper: 
	 *   
	 *   Automated Defect and Correlation Length Analysis of Block Copolymer Thin Film Nanopatterns
	 *   Jeffrey N. Murphy , Kenneth D. Harris, Jillian M. Buriak 
	 *   PLoS ONE 10(7): e0133088
	 *   Published: July 24, 2015
	 *   URL: https://doi.org/10.1371/journal.pone.0133088
	 *   DOI: 10.1371/journal.pone.0133088
	 */

// SETTINGS  /{{{
	program_name = "ADAblock";
	program_version = "v1.02"; prog_version = 1.02;
	modification_date = "2017.08.17";
	d_mode = 1; //(diagnostc mode)
	requires("1.49o"); // Requires Latest Version of ImageJ
	// http://fiji.sc/wiki/index.php/Auto_Threshold
	
	savestages = 0; /** Saves more images if set to 1 **/

// END OF SETTINGS  //}}}	


// MODIFIABLE DEFINED VARIABLES  //{{{
// Constants whose values can be modified.
image_log_name = "Output_Log";
period_limits_nm = newArray(10,150); // minimum & maximum period in nanometres
period_default_nm = 50;
auto_smoothing_factor = 0.15;
period_range_min_nm = 4; //nm // For FFT period auto-detection
period_range_max_nm = 100; //nm // For FFT period auto-detection
binary_grooming = true; // Include a binary grooming step to reduce extra defects.

// END OF MODIFIABLE DEFINED VARIABLES  //}}}

// DEPENDENCIES  //{{{
// Exit if necessary plugins are not installed.

// "Auto Local Threshold" Plugin
List.setCommands; 
if (List.get("Auto Local Threshold")!="") { 
	// plugin is installed 
} else {
	print("Install: AUTO LOCAL THRESHOLD");
	print("URL: http://bit.ly/plugin_ALT");
	exit("AUTO LOCAL THRESHOLD\nnot installed.");
}

// END OF DEPENDENCIES  //}}}

// FUNCTIONS  //{{{

// 000 Output Tags & Data //{{{
	var output_tags = newArray(0);
	var output_data= newArray(0);
	var output_labels = newArray(0);
	function outputTD(tag,data){
		output_tags = Array.concat(output_tags,tag); output_data = Array.concat(output_data,data);
	}
	function outputL(label){
		output_labels = Array.concat(output_labels,label);
	}
//}}}

// 001 AppendToArray  //{{{
// Appends a value to the array
// Returns the new array
// Format: array = AppendToArray(array,value);
function AppendToArray(array,value) {
	temp_array = newArray(lengthOf(array) + 1);
	for(i = 0; i < lengthOf(array); i++){
		temp_array[i] = array[i];
	}
	temp_array[lengthOf(temp_array) - 1] = value;
	array = temp_array;
	return array;
}  //}}}

// 002 CheckInArray  //{{{
// Checks if a value is already in the array
// Returns true or false
function CheckInArray(array,value) {
	i = 0; check = false;
	while(i < lengthOf(array)){
		if(array[i]==value){ i = lengthOf(array); check = true; }
		else{i++;}
	}
	return check;
}  //}}}

// 003 Indexer  //{{{
// Indexer: Converts (x,y) coordinates to values
// Only works for images up to 9999 x 9999
// Allows for creating a single list of values
function xyIndexer(x,y) {
	value = 100000000 + 10000*x + y;
	return value;
}  //}}}

// 004 DeIndexer  //{{{
// DeIndexer: Converts values to (x,y) coordinates
// Returns Array (x,y)
function xyDeIndexer(value) {
	x = floor((value - 100000000)/10000);
	y = value - 100000000 - x*10000;
	xy = newArray(x,y);
	return xy;
}  //}}}

// 005 Circuit Value  //{{{
// Circuit Value (Count the Jumps)
// Returns Number of Jumps Required (value)
// Note: outside of image, pixel value is 0
function circuitValue(x,y,value) {
	count = 0;
	if(getPixel(x+1,y+1)-getPixel(x,y+1)>=value){count++;}
	if(getPixel(x+1,y)-getPixel(x+1,y+1)>=value){count++;}
	if(getPixel(x+1,y-1)-getPixel(x+1,y)>=value){count++;}
	if(getPixel(x,y-1)-getPixel(x+1,y-1)>=value){count++;}
	if(getPixel(x-1,y-1)-getPixel(x,y-1)>=value){count++;}
	if(getPixel(x-1,y)-getPixel(x-1,y-1)>=value){count++;}
	if(getPixel(x-1,y+1)-getPixel(x-1,y)>=value){count++;}
	if(getPixel(x,y+1)-getPixel(x-1,y+1)>=value){count++;}
	return count;
}  //}}}

// 006 Four-Connected Circuit Value  //{{{
// Circuit Value (Count the Jumps), excluding diaconals (+)
// Returns Number of Jumps Required (only values: 0,1,2)
// Note: outside of image, pixel value is 0
function circuitValueFour(x,y,value) {
	count = 0;
	if(getPixel(x+1,y)-getPixel(x,y-1)>=value){count++;}
	if(getPixel(x,y-1)-getPixel(x-1,y)>=value){count++;}
	if(getPixel(x-1,y)-getPixel(x,y+1)>=value){count++;}
	if(getPixel(x,y+1)-getPixel(x+1,y)>=value){count++;}
	return count;
}

// 006b Four-Connected X
function circuitValueX(x,y,value) {
	count = 0;
	if(getPixel(x+1,y+1)-getPixel(x+1,y-1)>=value){count++;}
	if(getPixel(x+1,y-1)-getPixel(x-1,y-1)>=value){count++;}
	if(getPixel(x-1,y-1)-getPixel(x-1,y+1)>=value){count++;}
	if(getPixel(x-1,y+1)-getPixel(x+1,y+1)>=value){count++;}
	return count;
}  //}}}

// 007 Neighbour Value  //{{{
// Neighbour Value (or Point Value)
// Returns the number of 8-connected pixels (>= threshold value)
// Note: outside of image, pixel value is 0
function neighbourValue(x,y,value) {
	count = 0;
	if(getPixel(x,y+1)>=value){count++;}
	if(getPixel(x+1,y+1)>=value){count++;}
	if(getPixel(x+1,y)>=value){count++;}
	if(getPixel(x+1,y-1)>=value){count++;}
	if(getPixel(x,y-1)>=value){count++;}
	if(getPixel(x-1,y-1)>=value){count++;}
	if(getPixel(x-1,y)>=value){count++;}
	if(getPixel(x-1,y+1)>=value){count++;}
	return count;
}  //}}}

// 008 Neighbour Value Exact //{{{
// Neighbour Value (or Point Value)
// Returns the number of 8-connected pixels (>= threshold value)
// Note: outside of image, pixel value is 0
function neighbourValueExact(x,y,value) {
	count = 0;
	if(getPixel(x,y+1)==value){count++;}
	if(getPixel(x+1,y+1)==value){count++;}
	if(getPixel(x+1,y)==value){count++;}
	if(getPixel(x+1,y-1)==value){count++;}
	if(getPixel(x,y-1)==value){count++;}
	if(getPixel(x-1,y-1)==value){count++;}
	if(getPixel(x-1,y)==value){count++;}
	if(getPixel(x-1,y+1)==value){count++;}
	return count;
}

function neighbourValueExactFour(x,y,value) {
	count = 0;
	if(getPixel(x,y+1)==value){count++;}
	if(getPixel(x+1,y)==value){count++;}
	if(getPixel(x,y-1)==value){count++;}
	if(getPixel(x-1,y)==value){count++;}
	return count;
}  //}}}

// 009 FollowTwo  //{{{
// Follows a series of points in a skeletonized image
// Start at x,y Terminal Point; ends at Junction Point
// Requires ~4755 msec for "snake" of entire 1280x896 image
// Returns {x, y, distance travelled, pixels}
function followTwo(x,y,value) {
	condition = true; pixels = 0; u = x; v = y; a = 0; b = 0;
	while(condition){
		// Find which direction to go:
		count = 0; i = 0; j = 0; ip = u-x ; jp = v-y ; 
		// count = neighbours; i,j = direction to be moved; u,v = previous position; a = t-movement; b = x-movemnt;
		if(getPixel(x,y+1)>=value){count++; if(ip == 0 && jp == 1){ }else{ i = 0; j = 1; }}
		if(getPixel(x+1,y+1)>=value){count++; if(ip == 1 && jp == 1){ }else{ i = 1; j = 1; }}
		if(getPixel(x+1,y)>=value){count++; if(ip == 1 && jp == 0){ }else{ i = 1; j = 0; }}
		if(getPixel(x+1,y-1)>=value){count++; if(ip == 1 && jp == -1){ }else{ i = 1; j = -1; }}
		if(getPixel(x,y-1)>=value){count++; if(ip == 0 && jp == -1){ }else{ i = 0; j = -1; }}
		if(getPixel(x-1,y-1)>=value){count++; if(ip == -1 && jp == -1){ }else{ i = -1; j = -1; }}
		if(getPixel(x-1,y)>=value){count++; if(ip == -1 && jp == 0){ }else{ i = -1; j = 0; }}
		if(getPixel(x-1,y+1)>=value){count++; if(ip == -1 && jp == 1){ }else{ i = -1; j = 1; }}
		// setPixel(x,y,value); // print(steps + " " + x + " " + y);
		// Check to see if we should continue:
		if(pixels > 0 && count == 2){ u = x; v = y; x += i; y+= j; pixels++;} 
		else if(pixels==0 && count<=2){ u = x; v = y; x += i; y+= j; pixels++;}
		else {condition = false;  pixels++; i = 0; j = 0;}
		if(abs(i)+abs(j) == 1){a++;} else if(abs(i)+abs(j) == 2){b++;}
	}
	distance = a + sqrt(2)*b; result = newArray(x,y,distance,pixels);
	return result;
}  //}}}

// 010 FollowErase  //{{{
// Follows a series of points in a skeletonized image
// Start at x,y Terminal Point; deletes last unnecessary Junction Point
// Requires ~4708 msec for "snake" of entire 1280x896 image
// Returns {x, y, pixels}
function followErase(x,y,value) {
	condition = true; pixels = 0; erase = 10;
	while(condition){
		// Find which direction to go:
		ocount = 0; dcount = 0; i = 0; j = 0; 
		// count = neighbours; i,j = direction to be moved; u,v = previous position; a = t-movement; b = x-movemnt;
		if(getPixel(x,y+1)>=value){ocount++; i = 0; j = 1;}
		if(getPixel(x+1,y+1)>=value){dcount++; i = 1; j = 1;}
		if(getPixel(x+1,y)>=value){ocount++; i = 1; j = 0;}
		if(getPixel(x+1,y-1)>=value){dcount++; i = 1; j = -1;}
		if(getPixel(x,y-1)>=value){ocount++; i = 0; j = -1;}
		if(getPixel(x-1,y-1)>=value){dcount++; i = -1; j = -1;}
		if(getPixel(x-1,y)>=value){ocount++; i = -1; j = 0;}
		if(getPixel(x-1,y+1)>=value){dcount++; i = -1; j = 1;}
		count = ocount + dcount;
		if(count==1){setPixel(x,y,erase); pixels++; x += i; y+= j;}
		else {	condition = false;
			// Measure Circuit Value
			cvcount = 0; // Eight-connected
			if(getPixel(x+1,y+1)-getPixel(x,y+1)>=value){cvcount++;}
			if(getPixel(x+1,y)-getPixel(x+1,y+1)>=value){cvcount++;}
			if(getPixel(x+1,y-1)-getPixel(x+1,y)>=value){cvcount++;}
			if(getPixel(x,y-1)-getPixel(x+1,y-1)>=value){cvcount++;}
			if(getPixel(x-1,y-1)-getPixel(x,y-1)>=value){cvcount++;}
			if(getPixel(x-1,y)-getPixel(x-1,y-1)>=value){cvcount++;}
			if(getPixel(x-1,y+1)-getPixel(x-1,y)>=value){cvcount++;}
			if(getPixel(x,y+1)-getPixel(x-1,y+1)>=value){cvcount++;}
			// Check Conditions
			if(cvcount<2){setPixel(x,y,erase); pixels++;}
			else if(cvcount==2 && dcount==2){
				scount = 0; // Up,Left,Down,Right;
				if(getPixel(x,y+1)>=value){scount++; i = 0; j = 1;}
				if(getPixel(x+1,y)>=value){scount++; i = 1; j = 0;}
				if(getPixel(x,y-1)>=value){scount++; i = 0; j = -1;}
				if(getPixel(x-1,y)>=value){scount++; i = -1; j = 0;}
				setPixel(x,y,erase); pixels++;
			}
		}

	}
	result = newArray(x,y,pixels);
	return result;
}  //}}}

// 011 Set Foreground Index //{{{
// Sets the Foreground Colour using an Index
// Value corresponds to Color when using FloodFill
function setForegroundIndex(value) {
	getLut(reds,greens,blues);
	setForegroundColor(reds[value],greens[value],blues[value]);
} //}}}

// 012 LUT Functions //{{{

// LUTS One
function LUTs_001(true_or_false){
	if(true_or_false){
		reds = newArray(256); greens = newArray(256); blues = newArray(256);
		for(n=0; n<256; n++){ reds[n] = 255-n; greens[n] = 255-n; blues[n] = 255-n;}
		// Array.fill(reds, 0); Array.fill(greens, 0); Array.fill(blues, 0)
		value = 1; reds[value] = 255; greens[value] = 200; blues[value] = 200;
		value = 2; reds[value] = 120; greens[value] = 255; blues[value] = 120;
		value = 254; reds[value] = 0; greens[value] = 0; blues[value] = 170;
		value = 253; reds[value] = 0; greens[value] = 165; blues[value] = 120;
		setLut(reds, greens, blues);
	}
}


// LUT Index Modification 
function LUT_index_mod(index,r,g,b){
	getLut(reds, greens, blues);
	reds[index] = r; greens[index] = g; blues[index] = b;
	setLut(reds, greens, blues);
}


//}}}

// 013 Weighted Least Squares Algorithm //{{{
// For estimation of line widths

// Array Summing Function
function arraySumP3(ArrayA,ArrayB,ArrayC){
	sum = 0;
	if(ArrayB==1){ for(n=0; n<ArrayA.length; n++){ sum += ArrayA[n];} }
	else if(ArrayC==1){ for(n=0; n<ArrayA.length; n++){ sum += ArrayA[n]*ArrayB[n];} }
	else { for(n=0; n<ArrayA.length; n++){ sum += ArrayA[n]*ArrayB[n]*ArrayC[n];} }
	return sum;
}

// Simple Linear Regression (Weighted Least Squares)
// for measuring Average Width, using a series of Perimeters and Areas
// Formula: A = wP-C; A = area; P = perimeter; w = half-width; C = a constant
// Xi = Perimeters. Yi = Areas. (both are arrays)
// options: unweighted, inv_sqrt_area, median_distance
function particleWLSQ(Xi,Yi,Weighting_Method,iterations){
	N = Xi.length;

	// Weighting array
	Wi_not_norm = newArray(N); Wi = newArray(N);
	if(Weighting_Method=="unweighted"){ Array.fill(Wi_not_norm,1); }
	else if(Weighting_Method=="inv_sqrt_area"){ for(n=0; n<N; n++){ Wi_not_norm[n] = 1/(sqrt(Yi[n])); } }
	else if(Weighting_Method=="median_distance" && N >= 3){
		if(N%2==1){
			n_med = (N-1)/2+1; if(n_med>=Xi.length-1){n_med = Xi.length-2;}
			Xi_med = Xi[n_med+1]/3 + Xi[n_med]/3 + Xi[n_med-1]/3;
		}
		else{
			Xi_med = 0.5*Xi[floor(N/2)]+0.5*Xi[floor(N/2+1)];
		}
		for(n=0; n<N; n++){ 
			Wi_not_norm[n] = 1/pow(1+abs(Xi[n]-Xi_med),1.5); 
		}
	}
	else { Array.fill(Wi_not_norm,1);}
	Wi_Normalization = arraySumP3(Wi_not_norm,1,1);
	for(n=0; n<N; n++){ Wi[n] = Wi_not_norm[n] / Wi_Normalization; }
	Wi_method = Wi_not_norm;
	
	// Calculate All Sums
	S_Yi = arraySumP3(Yi,1,1);
	S_Xi = arraySumP3(Xi,1,1);
	S_Wi = arraySumP3(Wi,1,1);
	S_XiXi = arraySumP3(Xi,Xi,1);
	S_YiXi = arraySumP3(Yi,Xi,1);
	S_WiWi = arraySumP3(Wi,Wi,1);
	S_WiXi = arraySumP3(Wi,Xi,1);
	S_WiYi = arraySumP3(Wi,Yi,1);
	S_WiXiXi = arraySumP3(Wi,Xi,Xi);
	S_WiXiYi = arraySumP3(Wi,Xi,Yi);
	
	// Form of Y = aX+b
	Alpha = ( S_YiXi - (S_Yi*S_Xi/N)) / (S_XiXi - (S_Xi*S_Xi/N));
	Beta = (S_Yi/N) - Alpha*(S_Xi/N);
	AlphaW = (S_Wi*S_WiXiYi - S_WiXi*S_WiYi)/(S_Wi*S_WiXiXi - S_WiXi*S_WiXi);
	BetaW = (S_WiXiXi*S_WiYi - (S_WiXi*S_WiXiYi))/(S_Wi*S_WiXiXi - S_WiXi*S_WiXi);
	
	// Iterative Re-Weighting: Normalization
	BetaWW = BetaW; AlphaWW = AlphaW;
	for(m=0; m<iterations; m++){
		Ei = newArray(N); Y_fit = newArray(N);
		for(n=0; n<N; n++){
			Y_fit[n] = AlphaWW*Xi[n]+BetaWW;
			Ei[n] = abs(Y_fit[n] - Yi[n]);
			Wi_not_norm[n] = 1/Ei[n]*Wi_method[n];
		}
		Wi_Normalization = arraySumP3(Wi_not_norm,1,1);
		for(n=0; n<N; n++){ Wi[n] = Wi_not_norm[n] / Wi_Normalization; }
		
		// Re-Calculate Wi-Sums
		S_Wi = arraySumP3(Wi,1,1);
		S_WiWi = arraySumP3(Wi,Wi,1);
		S_WiXi = arraySumP3(Wi,Xi,1);
		S_WiYi = arraySumP3(Wi,Yi,1);
		S_WiXiXi = arraySumP3(Wi,Xi,Xi);
		S_WiXiYi = arraySumP3(Wi,Xi,Yi);
		
		// Form of Y = aX+b
		BetaWW = (S_WiXiXi*S_WiYi - (S_WiXi*S_WiXiYi))/(S_Wi*S_WiXiXi - S_WiXi*S_WiXi);
		AlphaWW = (S_Wi*S_WiXiYi - S_WiXi*S_WiYi)/(S_Wi*S_WiXiXi - S_WiXi*S_WiXi);	
	}
	W_Width_Calc = 2*AlphaWW;
	
	//Results: Alpha, Beta, AlphaW, BetaW, AlphaWW, BetaWW, W_Width_Calc
	results = newArray(W_Width_Calc, Alpha, Beta, AlphaW, BetaW, AlphaWW, BetaWW, iterations);
	return results;
}

// Array Maximum (not presently used)
// returns the maximum value of the array and the position in the array
function arrayMax(array){
	max_value = array[0]; i_max = 0;
	for(i=0; i<array.length; i++){
		if(array[i] > max_value){max_value = array[i]; i_max = i;}
	}
	result = newArray(i_max,max_value);
	return result;
}

// Array Maxima
// returns the maximum values of the array and the positions in the array
function arrayMaxima(Value_array,n_array,Value,n){
	for(i=0; i<Value_array.length; i++){
		if(Value_array[i]<Value){
			for(j=Value_array.length-1; j>i; j--){
				Value_array[j] = Value_array[j-1];
				n_array[j] = n_array[j-1];
			}
			Value_array[i] = Value;
			n_array[i] = n;
			result = Array.concat(Value_array,n_array);		
			return result;
		}
	}
	result = Array.concat(Value_array,n_array);		
	return result;
}

// Simple Linear Regression (Weighted Least Squares) **FULL**
// As above, except that Perimeters and Areas gathered from results table directly
// includes steps required for *exclusion* of small-area particles (min_area)  ## by area
// and exclusion via dropping the largest area values in the series (large_drop) ## by count
function particleWLSQfull(results_i,results_n,Weighting_Method,iterations,min_area,large_drops){
	count = -1*large_drops; drops = newArray(large_drops); Array.fill(drops, 0);
	n_drops = newArray(large_drops); Array.fill(n_drops, -1);
	if(large_drops>0){
		for(n=results_i; n<results_n; n++){ 
			area_value = getResult("Area",n);
			if(area_value >= min_area){
				count +=1;  
				if(area_value > drops[large_drops-1]){
					max_values = arrayMaxima(drops,n_drops,area_value,n);
					drops = Array.trim(max_values, large_drops);
					n_drops = Array.slice(max_values,large_drops,max_values.length);
				}
			}
		}
		if(count<2){ results = newArray(-1, -1, -1, -1, -1, -1, -1, -1); return results; }
		m = 0;
		Xi = newArray(count); Yi = newArray(count);
		for(n=results_i; n<results_n; n++){
			drops_onoff = true;
			for(nn=0; nn<drops.length; nn++){
				if(n_drops[nn]==n){drops_onoff = false;}
			}
			if((getResult("Area",n) >= min_area) & drops_onoff){
				Yi[m] = getResult("Area",n);
				Xi[m] = getResult("Perim.",n);
				m ++;
			}
		}		
	} else {
		for(n=results_i; n<results_n; n++){
			if(getResult("Area",n) >= min_area){
				count+= 1;
			}
		}
		if(count<3){ results = newArray(-1, -1, -1, -1, -1, -1, -1, -1); return results; }
		m = 0;
		Xi = newArray(count); Yi = newArray(count);
		for(n=results_i; n<results_n; n++){
			if(getResult("Area",n) >= min_area){
				Yi[m] = getResult("Area",n);
				Xi[m] = getResult("Perim.",n);
				m ++;
			}
		}
	}
	
	results = particleWLSQ(Xi,Yi,Weighting_Method,iterations);
	return results;
} //}}}

// 014 Edge Walk Pixels //{{{
// follows the edge of the image
// Checks to see if how many times the edge is touched by an object
// returns number of times touching and number of pixels touching

function edgeWalkPixels(colour_object,zero_x,zero_y,width,height) {
	jump = 1; j_count = 0; p_count = 0;
	px = getPixel(zero_x,zero_y); if(px==colour_object){old=1;}else{old=0;}
	top_count = 0; 
	y = zero_y;
	for(x=zero_x+1; x<width; x++){ px = getPixel(x,y); if(px==colour_object){new=1; p_count++;}else{new=0;} if(new-old==jump){j_count++;} old=new; }
	top_count = p_count;
	x = width-1; 
	for(y=zero_y+1; y<height; y++){ px = getPixel(x,y); if(px==colour_object){new=1; p_count++;}else{new=0;} if(new-old==jump){j_count++;} old=new; }
	right_count = p_count - top_count;
	y = height-1; 
	for(x=width-2; x>=zero_x; x--){ px = getPixel(x,y); if(px==colour_object){new=1; p_count++;}else{new=0;} if(new-old==jump){j_count++;} old=new; }
	bottom_count = p_count - (right_count + top_count);
	x = zero_x; 
	for(y=height-2; y>=zero_y; y--){ px = getPixel(x,y); if(px==colour_object){new=1; p_count++;}else{new=0;} if(new-old==jump){j_count++;} old=new; }
	left_count = p_count - (right_count + top_count + bottom_count);
	if(top_count>0){tc = 1000;}else{tc = 0;}
	if(right_count>0){rc = 100;}else{rc = 0;}
	if(bottom_count>0){bc = 10;}else{bc = 0;}
	if(left_count>0){lc = 1;}else{lc = 0;}
	sum = tc + rc + bc + lc; //!@#$ could used a more nuanced approach, but only an issue with small patterns.
	results = newArray(j_count,p_count,sum);
	return results;
}

// Edge Pixel Count: implements Edge Walk Pixels
// Requires results table. 
// No knowledge of phases is necessary
// except that colour_object cannot be the 
function edgePixelCount(colour_object,zero_x,zero_y,w,h){
	for(n=0; n<nResults; n++){
		xo = getResult("XStart",n); yo = getResult("YStart",n); xy_value = getPixel(xo,yo);
		bx = getResult("BX",n); by = getResult("BY",n);
		Lx = getResult("Width",n); Ly = getResult("Height",n);
		if(bx==0 || by==0 || bx+Lx==w || by+Ly==h){
			setResult("OnEdge",n,1);
			setForegroundIndex(colour_object);
			floodFill(xo,yo);
			result = edgeWalkPixels(colour_object,zero_x,zero_y,w,h);
			setResult("EdgeTouch",n,result[0]);
			setResult("EdgePixels",n,result[1]);
			setResult("Sides",n,result[2]);
			//if(n<nPositive){setForegroundIndex(255);}else{setForegroundIndex(0);}
			setForegroundIndex(xy_value); // Replacement line
			floodFill(xo,yo);
		}
		else{ setResult("OnEdge",n,0); }
	}
	updateResults(); return 1;
} //}}}

// 015 CHECK INSIDE //{{{
// Discover whether objects are enclosed inside of another. Enclosure indicates exterior lines or a defect
// contains particles: "Contains" = 1
// enclosed by particles: "Enclosed"

function checkInside(nPositive,nTotal,cvalue){
	if(isNaN(getResult("Contains",0))){for(n=0; n<nTotal;n++){setResult("Contains",n,0);}}
	if(isNaN(getResult("Enclosed",0))){for(n=0; n<nTotal;n++){setResult("Enclosed",n,0);}}
	if(isNaN(getResult("Enclosed.By",0))){for(n=0; n<nTotal;n++){setResult("Enclosed.By",n,-1);}}
	
	for(i=0; i<2; i++){
		if(i==0){start_A = 0; end_B = nPositive; start_C = nPositive; end_D = nTotal;}
		else{start_A = nPositive; end_B = nTotal; start_C = 0; end_D = nPositive;}
		for(n=start_A; n<end_B; n++){
			m_inside_check = newArray(0);
			for(m=start_C; m<end_D; m++){
				nBx = getResult("BX",n); nBy = getResult("BY",n);
				nLx = getResult("Width",n); nLy = getResult("Height",n);
				mBx = getResult("BX",m); mBy = getResult("BY",m);
				mLx = getResult("Width",m); mLy = getResult("Height",m);
				// if "m" entry is inside of "n" entry
				if(mBx>nBx && mBy>nBy && (mBx+mLx)<(nBx+nLx) && (mBy+mLy)<(nBy+nLy)){
					// Checks to see if neighbouring values changed after flood fill. Indicates enclosure.
					// Enclosed ones, via ">" (not ">=" definition have first pixel inside regardless.
					m_inside_check = Array.concat(m_inside_check,m);
				}
				//else{if(getResult("Lines",n)!=1){setResult("Lines",n,0);}}
			}
			// Flood Fill part is moved outside for speed...
			if(m_inside_check.length>0){
				nx = getResult("XStart",n); ny = getResult("YStart",n); nxy_value = getPixel(nx,ny);
				setForegroundIndex(cvalue); floodFill(nx,ny,"8-connected");
				for(m=0; m<m_inside_check.length; m++){
					mx = getResult("XStart",m_inside_check[m]); my = getResult("YStart",m_inside_check[m]);
					mn_connection = neighbourValueExact(mx,my,cvalue);
					if(mn_connection>0){
						count = getResult("Contains",n); setResult("Contains",n,count+1);
						setResult("Enclosed",m_inside_check[m],1); setResult("Enclosed.By",m_inside_check[m],n);
					}
				}
				setForegroundIndex(nxy_value); floodFill(nx,ny,"8-connected"); // return to original value
			}
		}
	}
	updateResults(); return 1;
} 

function checkInsideEDGE(nPositive,nTotal,cvalue,w,h){
	if(isNaN(getResult("Contains",0))){for(n=0; n<nTotal;n++){setResult("Contains",n,0);}}
	if(isNaN(getResult("Enclosed",0))){for(n=0; n<nTotal;n++){setResult("Enclosed",n,0);}}
	if(isNaN(getResult("Enclosed.By",0))){for(n=0; n<nTotal;n++){setResult("Enclosed.By",n,-1);}}
	
	for(i=0; i<2; i++){
		if(i==0){start_A = 0; end_B = nPositive; start_C = nPositive; end_D = nTotal;}
		else{start_A = nPositive; end_B = nTotal; start_C = 0; end_D = nPositive;}
		for(n=start_A; n<end_B; n++){
			m_inside_check = newArray(0);
			for(m=start_C; m<end_D; m++){
				nBx = getResult("BX",n); nBy = getResult("BY",n);
				nLx = getResult("Width",n); nLy = getResult("Height",n);
				mBx = getResult("BX",m); mBy = getResult("BY",m);   
				mLx = getResult("Width",m); mLy = getResult("Height",m);
				// if "m" entry is inside of "n" entry
				if(mBx>=nBx && mBy>=nBy && (mBx+mLx)<=(nBx+nLx) && (mBy+mLy)<=(nBy+nLy)){
					// Checks to see if neighbouring values changed after flood fill. Indicates enclosure.
					// Enclosed ones, via ">" (not ">=" definition have first pixel inside regardless.
					m_inside_check = Array.concat(m_inside_check,m);
				}
				//else{if(getResult("Lines",n)!=1){setResult("Lines",n,0);}}
			}
			// Flood Fill part is moved outside for speed...
			if(m_inside_check.length>0){
				nx = getResult("XStart",n); ny = getResult("YStart",n); nxy_value = getPixel(nx,ny);
				setForegroundIndex(cvalue); floodFill(nx,ny,"8-connected");
				if( (nx==0 || nx==w-1) || (ny==0 || ny==h-1) ){ nxny = furthestNonEdgePixel(nx,ny,w,h,cvalue); nx = nxny[0]; ny = nxny[1];}
				for(m=0; m<m_inside_check.length; m++){
					mx = getResult("XStart",m_inside_check[m]); my = getResult("YStart",m_inside_check[m]);
					mn_connection = neighbourValueExact(mx,my,cvalue);
					if(mn_connection>0){
						count = getResult("Contains",n); setResult("Contains",n,count+1);
						setResult("Enclosed",m_inside_check[m],1); setResult("Enclosed.By",m_inside_check[m],n);
					}
				}
				setForegroundIndex(nxy_value); floodFill(nx,ny,"8-connected"); // return to original value
			}
		}
	}
	updateResults(); return 1;
}//}}}

// 016 Conditions True //{{{
// If x conditions are *true*, then 
// e.g. 
// array_labels = newArray("Area","Circ.");
// array_conditions = newArray(">=",">");
// array_values = newArray(50,0.80);
// column is string to name COLUMN in Results Table
// n_start = 0; n_end = nResults

function conditionsTrue(array_labels,array_conditions,array_values,n_start,n_end,column,x){
	// array_conditions: ==,!=,>,<,>=,<=
	product_sum = 0; sum_sum = 0;
	cond_eval = newArray(array_values.length);
	for(n=n_start; n<n_end; n++){
		for(m=0; m<array_values.length; m++){
			label = array_labels[m];
			cond = array_conditions[m];
			value = array_values[m];
			L_value = getResult(label,n); 
			macro_expression = "result="+L_value+cond+value+"; return toString(result);";
			cond_eval[m] = eval(macro_expression);
		}
		product = 1; sum = 0;
		for(m=0; m<cond_eval.length; m++){
			product = product * cond_eval[m];
			sum += cond_eval[m];
		}
		if(sum >= x){result = 1; sum_sum += 1;}else{result = 0;}
		setResult(column,n,result);
		product_sum += product;
	}
	updateResults(); return sum_sum;
}

// Version with weighs attached
// prod_true = array of values if true
// prod_false = array of values if true
// binary = 1 or 0. if true, results will be 1 or 0; if false, they will be 
function conditionsWeighted(array_labels,array_conditions,array_values,prod_true,prod_false,binary,threshold,n_start,n_end,column){
	// array_conditions: ==,!=,>,<,>=,<=
	product_sum = 0; sum_sum = 0;
	cond_eval = newArray(array_values.length);
	for(n=n_start; n<n_end; n++){
		for(m=0; m<array_values.length; m++){
			label = array_labels[m];
			cond = array_conditions[m];
			value = array_values[m];
			L_value = getResult(label,n); 
			macro_expression = "result="+L_value+cond+value+"; return toString(result);";
			if(eval(macro_expression)){cond_eval[m] = prod_true[m]}else{cond_eval[m] = prod_false[m];}
		}
		product = 1; sum = 0;
		for(m=0; m<cond_eval.length; m++){
			product = product * cond_eval[m];
			sum += cond_eval[m];
		}
		if(binary){
			if(sum >= threshold){result = 1; sum_sum += 1;}else{result = 0;}
		} else {
			if(sum >= threshold){sum_sum += 1;}
			result = sum;
		}
			
		setResult(column,n,result);
		product_sum += product;
	}
	updateResults(); return sum_sum;
} 

example = 0;
if(example){
	start = getTime();
	array_labels = newArray("Phase","Area","Circ.");
	array_conditions = newArray("==",">=",">");
	array_values = newArray(1,50,0.72);
	prod_true = newArray(1,1,1);
	prod_false = newArray(0,0,0);
	n_start = 0;
	n_end = nResults();
	x = 3;
	threshold = 3;
	binary = 0;
	columnA = "CT";
	columnB = "CW";
	
	a = conditionsTrue(array_labels,array_conditions,array_values,n_start,n_end,columnA,x);
	print("A: " + a);
	b =  conditionsWeighted(array_labels,array_conditions,array_values,prod_true,prod_false,binary,threshold,n_start,n_end,columnB);
	print("B: " + b);
	end = getTime();
	print("Time: " + end-start);
} //}}}

// 017 Colour Particles //{{{
// Two functions for colouring particles according to a SINGLE condition.
// latter function includes PHASE SELECTION

// e.g.
//cP = colourParticles("Phase","==",1,100,0,nResults);
//cP = colourParticlesPhase(1,"Area","<",100,100,0,nResults);

function colourParticles(column,condition,value,index,n_start,n_end){
	count = 0; setForegroundIndex(index);
	for(n=n_start; n<n_end; n++){
		col_value = getResult(column,n);
		macro_expression = "result="+col_value+condition+value+"; return toString(result);";
		if(eval(macro_expression)){
			x = getResult("XStart",n);
			y = getResult("YStart",n);
			floodFill(x,y,"8-connected");
			count++;
		}	
	}
	return count;
}



function colourParticlesPhase(phase,column,condition,value,index,n_start,n_end){
	count = 0; setForegroundIndex(index);
	for(n=n_start; n<n_end; n++){
		if(getResult("Phase",n)==phase){
			col_value = getResult(column,n);
			macro_expression = "result="+col_value+condition+value+"; return toString(result);";
			if(eval(macro_expression)){
				x = getResult("XStart",n);
				y = getResult("YStart",n);
				floodFill(x,y,"8-connected");
				count++;
			}
		}
	}
	return count;
} 
/* INCOMPLETE
function colourByParameter(column,value,index,n_start,n_end){
	count = 0; setForegroundIndex(index);
	for(n=n_start; n<n_end; n++){
		col_value = getResult(column,n);
		macro_expression = "result="+col_value+condition+value+"; return toString(result);";
		if(eval(macro_expression)){
			x = getResult("XStart",n);
			y = getResult("YStart",n);
			floodFill(x,y,"8-connected");
			count++;
		}	
	}
	return count;
}
 */

//}}}

// 018 XY CODER (Encoder-Decoder) //{{{
/** 	type = "enc" will encode x & y to return value (just leave value = 0)
	type = "dec" will decode value to return x & y (just leave x & y = 0) **/
function xy_coder(x,y,value,type){
	scale = 10000;
	if(type=="enc"){
		enc = scale*x+y;
		return enc;
	}
	else if(type=="dec"){
		x = floor(value/scale);
		y = value-x*scale;
		dec = newArray(x,y);
		return dec;
	}
}


//}}}

// 019 Centre Pixel //{{{
/** 	Takes an array of pixels & returns centre-most pixel
	enc = 1 will encode; enc = 0 will give array **/
function centrePixel(xpoints,ypoints,enc) {
	xsum = 0; ysum = 0;
	for(i=0; i<xpoints.length; i++){ 
		xsum += xpoints[i]; 
		ysum += ypoints[i]; 
	}
	xavg = ysum / xpoints.length; yavg = ysum / ypoints.length;
	xo = xpoints[0]; yo = ypoints[0];
	min_dist = sqrt(pow(xavg-xo,2)+pow(yavg-yo,2));
	for(i=1; i<xpoints.length; i++){
		dist = sqrt(pow(xavg-xpoints[i],2)+pow(yavg-ypoints[i],2));
		if(dist<min_dist){
			xo = xpoints[i];
			yo = ypoints[i];
			min_dist = dist;
		}
	}
	if(enc == 1){ xy = xy_coder(xo,yo,0,"enc"); return xy;}
	else{ xy = newArray(xo,yo); return xy;	}
}
//}}}

// 020 IsPointInPath //{{{
/** 	Adapted from Python code for the "EVEN-ODD RULE"
	source: http://en.wikipedia.org/wiki/Even-odd_rule  **/


function isPointInPath(x, y, xpoly, ypoly){
	num = xpoly.length;
	i = 0;
	j = num - 1;
	c = false;
	for(i=0; i<num; i++){
		if( ((ypoly[i] > y) != (ypoly[j] > y)) & (x < (xpoly[j] - xpoly[i]) * (y - ypoly[i]) / (ypoly[j] - ypoly[i]) + xpoly[i]) ){
			if(c){ c = false; } else { c = true; }  // c = not c
		}
		j = i;
	}
	return c;
}



// Example
// Make selection of an object prior to running this
// Object should fill with set value.
/**
getSelectionCoordinates(xpoints,ypoints);
getSelectionBounds(xo, yo, wS, hS);  // wS & hS are width & height of selection
for(y = yo; y < yo+hS; y++){
	for(x = xo; x < xo+wS; x++){
		a = isPointInPath(x,y,xpoints,ypoints);
		if(a){setPixel(x,y,150);}
	}
}
**/

/** Original PYTHON CODE:

def isPointInPath(x, y, poly):
    num = len(poly)
    i = 0
    j = num - 1
    c = False
    for i in range(num):
        if  ((poly[i][1] > y) != (poly[j][1] > y)) and (x < (poly[j][0] - poly[i][0]) * (y - poly[i][1]) / (poly[j][1] - poly[i][1]) + poly[i][0]):
            c = not c
        j = i
    return c

**/

//}}}

// 021 Get Selection Pixels //{{{
// getSelectionPixels
/** 	xORy = 0 gives xpoints
	xORy = 1 gives ypoints
	xORy = 2 gives encoded points (xy_coder)
	xORy = 3 gives centre pixel, encoded
	xORy = 4 gives centre pixel, array **/
function getSelectionPixels(xORy){
	getSelectionCoordinates(xpoly,ypoly);
	getSelectionBounds(xo, yo, wS, hS);
	i = 0;
	for(y = yo; y < yo+hS; y++){
		for(x = xo; x < xo+wS; x++){
			if(isPointInPath(x,y,xpoly,ypoly)){i++;}
		}
	}
	xpoints = newArray(i); ypoints = newArray(i); i = 0;
	for(y = yo; y < yo+hS; y++){
		for(x = xo; x < xo+wS; x++){
			if(isPointInPath(x,y,xpoly,ypoly)){xpoints[i] = x; ypoints[i] = y; i++;}
		}
	}
	if(xORy==0){ return xpoints; }
	if(xORy==1){ return ypoints; }
	if(xORy==2){
		xy_points = newArray(xpoints.length);
		for(i=0; i<xpoints.length; i++){
			xy_points[i] = xy_coder(xpoints[i],ypoints[i],0,"enc");
		}
		return xy_points;
	}
	if(xORy==3){
		xy = centrePixel(xpoints,ypoints,1);
		return xy;
	}
	if(xORy==4){
		xy = centrePixel(xpoints,ypoints,0);
		return xy;
	}		
}
//}}}

// 022 Defect Encoder / Decoder  //{{{
/**	type = "enc" : encode. Input = phase, connectivity, x, y.
	type = "dec" : decode. Input = values. **/
function defect_coder(values, phase, connectivity, x, y,type){
	// encode values: 1PCCXXXXXYYYYY 
	// Phase (1,0); Connectivity (0...9); X-coordinates; Y-coordinates
	expts = newArray(12,10,5,0);
	if(type=="enc"){
		values = newArray(phase, connectivity, x, y);
		enc = 1 * pow(10,expts[0]+1);
		for(i=0; i<4; i++){ enc += values[i] * pow(10,expts[i]); }
		return enc;
	}
	else if(type=="dec"){
		dec = newArray(4);
		v = values;
		vs = 1 * pow(10,expts[0]+1);
		for(i=0; i<4; i++){
			v = v - vs;
			dec[i] = floor(v / pow(10,expts[i]));
			vs = dec[i] * pow(10,expts[i]);
		}
		return dec;
	}
}  //}}}

// 023 Integer String //{{{
/** converts an integer to a string for printing **/
function integerString(n){
	m = floor(log(n)/log(10));
	str = "";
	rem = n;
	for(i=m; i>=0; i--){
		ni = floor(rem/pow(10,i));
		rem = rem - ni * pow(10,i);
		str += toString(ni);
	}
	return str;
}
//}}}

// 024 Shift Values //{{{
/**	how == 0 use xpts & ypts pixel arrays
	how = 1 obtain xpts & ypts from selection 
	shift = how much the values will be shifted by
	e.g. 100: values will be shifted by +100
	-50: values will be shifted by -50	**/
function shiftValues(xpoints,ypoints,shift,how){
	if(how==1){
		getSelectionCoordinates(xpoly,ypoly);
		getSelectionBounds(xo, yo, wS, hS);
		i = 0;
		for(y = yo; y < yo+hS; y++){
			for(x = xo; x < xo+wS; x++){
				if(isPointInPath(x,y,xpoly,ypoly)){i++;}
			}
		}
		xpoints = newArray(i); ypoints = newArray(i); i = 0;
		for(y = yo; y < yo+hS; y++){
			for(x = xo; x < xo+wS; x++){
				if(isPointInPath(x,y,xpoly,ypoly)){xpoints[i] = x; ypoints[i] = y; i++;}
			}
		}
	}
	for(i=0; i<xpoints.length; i++){
		x = xpoints[i]; y = ypoints[i];
		value = getPixel(x,y);
		setPixel(x,y,value+shift);
	}
	return xpoints.length;
} //}}}

// 025 Distance from Edge //{{{
function edgeDistance(w,h,x,y){
	dist = minOf(minOf(x,w-x-1),minOf(y,h-y-1));
	return dist;
}//}}}


// 025 Junction Degree //{{{
// Junction degree search
/** 	value = pixel value sought
	newvalue = replacement pixel value	**/
function junctionDegree(value,newvalue){
	getSelectionCoordinates(xpoly,ypoly);
	getSelectionBounds(xo, yo, wS, hS);
	i = 0;
	for(y = yo; y < yo+hS; y++){
		for(x = xo; x < xo+wS; x++){
			if(isPointInPath(x,y,xpoly,ypoly)){i++;}
		}
	}
	xpoints = newArray(i); ypoints = newArray(i); i = 0;
	for(y = yo; y < yo+hS; y++){
		for(x = xo; x < xo+wS; x++){
			if(isPointInPath(x,y,xpoly,ypoly)){xpoints[i] = x; ypoints[i] = y; i++;}
		}
	}
	count = 0;
	for(i=0; i<xpoints.length; i++){
		x = xpoints[i]; y = ypoints[i];
		if(getPixel(x,y+1)==value){count++; setPixel(x,y+1,newvalue);}
		if(getPixel(x+1,y+1)==value){count++; setPixel(x+1,y+1,newvalue);}
		if(getPixel(x+1,y)==value){count++; setPixel(x+1,y,newvalue);}
		if(getPixel(x+1,y-1)==value){count++; setPixel(x+1,y-1,newvalue);}
		if(getPixel(x,y-1)==value){count++; setPixel(x,y-1,newvalue);}
		if(getPixel(x-1,y-1)==value){count++; setPixel(x-1,y-1,newvalue);}
		if(getPixel(x-1,y)==value){count++; setPixel(x-1,y,newvalue);}
		if(getPixel(x-1,y+1)==value){count++; setPixel(x-1,y+1,newvalue);}
	}
	return count;	
} //}}}

// 026 LER Functions //{{{
/** from LER_new_multi_lines_06_20140124.ijm **/
  ////////////////////////////////////////////
  // FUNCTIONS                      
  ////////////////////////////////////////////
  
  // FUNCTIONS
  // Point Value
  function PTV(x,y) {
    ptvalue = getPixel(x-1,y-1)+getPixel(x-1,y)+getPixel(x-1,y+1)+getPixel(x,y+1)+getPixel(x+1,y+1)+getPixel(x+1,y)+getPixel(x+1,y-1)+getPixel(x,y-1);
    ptv = ptvalue/255;
    return ptv;
  }
  
  // Normalized Dot Product: a.b/(|a||b|)
  function NDP(vecxa,vecya,vecxb,vecyb){
    A = sqrt(vecxa*vecxa+vecya*vecya);
    B = sqrt(vecxb*vecxb+vecyb*vecyb);
    DP = vecxa*vecxb+vecya*vecyb;
    return DP/(A*B);
  }
  
  // OrthoProject Function (tells you the orthogonal vector required to complete the triangle).
  var c_vector = newArray(2);
  // a & b must be arrays; xy-vectors
  function OrthoProject(a,b) {
    La = sqrt(pow(a[0],2)+pow(a[1],2));
    Lb = sqrt(pow(b[0],2)+pow(b[1],2));
    // unit_b = newArray(2); unit_b[0] = b[0]/Lb; unit_b[1] = b1]/Lb;
    value = (a[0]*b[0]+a[1]*b[1])/(La*Lb); if(value>1){value=1;}
    theta=acos(value);
    //print(theta*180/PI);
    s = Lb*value; f = La/s;
    c_vector = newArray(2); c_vector[0] = (f*b[0]-a[0])/2; c_vector[1] = (f*b[1]-a[1])/2;
    Lc = sqrt(pow(c_vector[0],2)+pow(c_vector[1],2));
    return Lc;
  }
  
  // Cross Product
  function CrossProdn(a,b) {
    zc = a[0]*b[1]-a[1]*b[0];
    if(zc==0){zcn=0;}else{zcn = zc/abs(zc);}
    return zcn;
  }
  
  // Point Value Function (Reverse Method = counting blanks)
  function PTVr(x,y) {
        if(getPixel(x-1,y-1)==0){p1=1;}else{p1=0;}
        if(getPixel(x-1,y)==0){p2=1;}else{p2=0;}
        if(getPixel(x-1,y+1)==0){p3=1;}else{p3=0;}
        if(getPixel(x,y+1)==0){p4=1;}else{p4=0;}
        if(getPixel(x+1,y+1)==0){p5=1;}else{p5=0;}
        if(getPixel(x+1,y)==0){p6=1;}else{p6=0;}
        if(getPixel(x+1,y-1)==0){p7=1;}else{p7=0;}
        if(getPixel(x,y-1)==0){p8=1;}else{p8=0;}
    ptv_blank = p1+p2+p3+p4+p5+p6+p7+p8;
    ptvr = 8-ptv_blank;
    return ptvr;
  }
  
  function PTVrPlus(x,y) {
        //if(getPixel(x-1,y-1)==0){p1=1;}else{p1=0;}
        if(getPixel(x-1,y)==0){p2=1;}else{p2=0;}
        //if(getPixel(x-1,y+1)==0){p3=1;}else{p3=0;}
        if(getPixel(x,y+1)==0){p4=1;}else{p4=0;}
        //if(getPixel(x+1,y+1)==0){p5=1;}else{p5=0;}
        if(getPixel(x+1,y)==0){p6=1;}else{p6=0;}
        //if(getPixel(x+1,y-1)==0){p7=1;}else{p7=0;}
        if(getPixel(x,y-1)==0){p8=1;}else{p8=0;}
    ptv_blank = p2+p4+p6+p8;
    ptvr = 4-ptv_blank;
    return ptvr;
  }
  
  function PTVmp(x,y,MaxPx) {
        if(getPixel(x-1,y-1)==MaxPx){p1=1;}else{p1=0;}
        if(getPixel(x-1,y)==MaxPx){p2=1;}else{p2=0;}
        if(getPixel(x-1,y+1)==MaxPx){p3=1;}else{p3=0;}
        if(getPixel(x,y+1)==MaxPx){p4=1;}else{p4=0;}
        if(getPixel(x+1,y+1)==MaxPx){p5=1;}else{p5=0;}
        if(getPixel(x+1,y)==MaxPx){p6=1;}else{p6=0;}
        if(getPixel(x+1,y-1)==MaxPx){p7=1;}else{p7=0;}
        if(getPixel(x,y-1)==MaxPx){p8=1;}else{p8=0;}
    ptv = p1+p2+p3+p4+p5+p6+p7+p8;
    return ptv;
  }
  
  function Follow(x,y,xo,yo,MaxPx){
    ip = xo-x; jp = yo-y;
    if(getPixel(x-1,y-1)==MaxPx){ if(ip==-1&&jp==-1){}else{walk=0;}}
    if(getPixel(x,y-1)==MaxPx){ if(ip==0&&jp==-1){}else{walk=1;}}
    if(getPixel(x+1,y-1)==MaxPx){ if(ip==1&&jp==-1){}else{walk=2;}}
    if(getPixel(x+1,y)==MaxPx){ if(ip==1&&jp==0){}else{walk=3;}}
    if(getPixel(x+1,y+1)==MaxPx){ if(ip==1&&jp==1){}else{walk=4;}}
    if(getPixel(x,y+1)==MaxPx){ if(ip==0&&jp==1){}else{walk=5;}}
    if(getPixel(x-1,y+1)==MaxPx){ if(ip==-1&&jp==1){}else{walk=6;}}
    if(getPixel(x-1,y)==MaxPx){ if(ip==-1&&jp==0){}else{walk=7;}}
    return walk;
  }
    
    i_walk = newArray(-1,0,1,1,1,0,-1,-1);
    j_walk = newArray(-1,-1,-1,0,1,1,1,0);
  
    // Phase
  function Phase(x,y,MaxPx,MinPx){
    Up = 0; Down = 0;
    if(getPixel(x-1,y-1)==MaxPx){a=1;} else if(getPixel(x-1,y-1)==MinPx){a=0;} 
    if(getPixel(x,y-1)==MaxPx){b=1;} else if(getPixel(x,y-1)==MinPx){b=0;}
    if(getPixel(x+1,y-1)==MaxPx){c=1;} else if(getPixel(x+1,y-1)==MinPx){c=0;} 
    if(getPixel(x+1,y)==MaxPx){d=1;} else if(getPixel(x+1,y)==MinPx){d=0;}
    if(getPixel(x+1,y+1)==MaxPx){e=1;} else if(getPixel(x+1,y+1)==MinPx){e=0;} 
    if(getPixel(x,y+1)==MaxPx){f=1;} else if(getPixel(x,y+1)==MinPx){f=0;}
    if(getPixel(x-1,y+1)==MaxPx){g=1;} else if(getPixel(x-1,y+1)==MinPx){g=0;}
    if(getPixel(x-1,y)==MaxPx){h=1;} else if(getPixel(x-1,y)==MinPx){h=0;} 
    if(b-a==1){Up++;}else if(b-a==-1){Down++;}
    if(c-b==1){Up++;}else if(c-b==-1){Down++;}
    if(d-c==1){Up++;}else if(d-c==-1){Down++;}
    if(e-d==1){Up++;}else if(e-d==-1){Down++;}
    if(f-e==1){Up++;}else if(f-e==-1){Down++;}
    if(g-f==1){Up++;}else if(g-f==-1){Down++;}
    if(h-g==1){Up++;}else if(h-g==-1){Down++;}
    if(a-h==1){Up++;}else if(a-h==-1){Down++;}
    return Up;
  }

  // NEW FUNCTIONS






function followTwoDistances(x,y,value) {
	condition = true; pixels = 0; u = x; v = y; a = 0; b = 0;
	distance_array = newArray(0);
	while(condition){
		distance = a + sqrt(2)*b;
		distance_array = Array.concat(distance_array,distance);
		// Find which direction to go:
		count = 0; i = 0; j = 0; ip = u-x ; jp = v-y ;
		// count = neighbours; i,j = direction to be moved; u,v = previous position; a = t-movement; b = x-movemnt;
		if(getPixel(x,y+1)>=value){count++; if(ip == 0 && jp == 1){ }else{ i = 0; j = 1; }}
		if(getPixel(x+1,y+1)>=value){count++; if(ip == 1 && jp == 1){ }else{ i = 1; j = 1; }}
		if(getPixel(x+1,y)>=value){count++; if(ip == 1 && jp == 0){ }else{ i = 1; j = 0; }}
		if(getPixel(x+1,y-1)>=value){count++; if(ip == 1 && jp == -1){ }else{ i = 1; j = -1; }}
		if(getPixel(x,y-1)>=value){count++; if(ip == 0 && jp == -1){ }else{ i = 0; j = -1; }}
		if(getPixel(x-1,y-1)>=value){count++; if(ip == -1 && jp == -1){ }else{ i = -1; j = -1; }}
		if(getPixel(x-1,y)>=value){count++; if(ip == -1 && jp == 0){ }else{ i = -1; j = 0; }}
		if(getPixel(x-1,y+1)>=value){count++; if(ip == -1 && jp == 1){ }else{ i = -1; j = 1; }}
		// setPixel(x,y,value); // print(steps + " " + x + " " + y);
		// Check to see if we should continue:
		if(pixels > 0 && count == 2){ u = x; v = y; x += i; y+= j; pixels++;} 
		else if(pixels==0 && count<=2){ u = x; v = y; x += i; y+= j; pixels++;}
		else {condition = false;  pixels++; i = 0; j = 0; }
		if(abs(i)+abs(j) == 1){a++;} else if(abs(i)+abs(j) == 2){b++;}

	}
	distance = a + sqrt(2)*b; result = newArray(x,y,distance,pixels);
	return distance_array;
}



function followTwoEncode(x,y,value) {
	condition = true; pixels = 0; u = x; v = y; a = 0; b = 0;
	encoded = newArray(0);
	while(condition){
		// Find which direction to go:
		count = 0; i = 0; j = 0; ip = u-x ; jp = v-y ; 
		// count = neighbours; i,j = direction to be moved; u,v = previous position; a = t-movement; b = x-movemnt;
		if(getPixel(x,y+1)>=value){count++; if(ip == 0 && jp == 1){ }else{ i = 0; j = 1; }}
		if(getPixel(x+1,y+1)>=value){count++; if(ip == 1 && jp == 1){ }else{ i = 1; j = 1; }}
		if(getPixel(x+1,y)>=value){count++; if(ip == 1 && jp == 0){ }else{ i = 1; j = 0; }}
		if(getPixel(x+1,y-1)>=value){count++; if(ip == 1 && jp == -1){ }else{ i = 1; j = -1; }}
		if(getPixel(x,y-1)>=value){count++; if(ip == 0 && jp == -1){ }else{ i = 0; j = -1; }}
		if(getPixel(x-1,y-1)>=value){count++; if(ip == -1 && jp == -1){ }else{ i = -1; j = -1; }}
		if(getPixel(x-1,y)>=value){count++; if(ip == -1 && jp == 0){ }else{ i = -1; j = 0; }}
		if(getPixel(x-1,y+1)>=value){count++; if(ip == -1 && jp == 1){ }else{ i = -1; j = 1; }}
		// setPixel(x,y,value); // print(steps + " " + x + " " + y);
		// Check to see if we should continue:
		encoded = Array.concat(encoded,xy_coder(x,y,0,"enc"));
		if(pixels > 0 && count == 2){ u = x; v = y; x += i; y+= j; pixels++;} 
		else if(pixels==0 && count<=2){ u = x; v = y; x += i; y+= j; pixels++;}
		else {condition = false;  pixels++; i = 0; j = 0; }
		if(abs(i)+abs(j) == 1){a++;} else if(abs(i)+abs(j) == 2){b++;}

	}
	distance = a + sqrt(2)*b; result = newArray(x,y,distance,pixels);
	return encoded;
}




// Note: outside of image, pixel value is 0
  
function findFirst(x,y,minvalue,findvalue) {
	u_array = newArray(-1,-1,-1,0,1,1,1,0);
	v_array = newArray(-1,0,1,1,1,0,-1,-1);
	x_list = newArray(1); x_list[0] = x;
	y_list = newArray(1); y_list[0] = y;
	i = 0;
	while(i<x_list.length){
		x = x_list[i]; y = y_list[i]; j = 0;
		if(getPixel(x,y)==findvalue){return newArray(x,y);}
		while(j<8){
			u = x + u_array[j]; v = y + v_array[j]; //print(u + " " + v);
			if(getPixel(u,v)>=minvalue){
				if(CheckInArrayXY(x_list,y_list,u,v)){ }
				else{
					x_list = Array.concat(x_list,u); y_list = Array.concat(y_list,v); 
					if(getPixel(u,v)==findvalue){return newArray(u,v);}
				}
			}
			j++;
		}
		oldvalue = getPixel(x,y);
		//setPixel(x,y,findvalue-1);
		//print(i + " " + oldvalue);
		i++;
	}
	return newArray(-1,-1);
}  

// CheckInArrayXY  //{{{
// Returns true or false
function CheckInArrayXY(x_array,y_array,x_value,y_value) {
	i = 0; check = false;
	while(i < lengthOf(x_array)){
		if(x_array[i]==x_value){
			if(y_array[i]==y_value){
				i = lengthOf(x_array); check = true;
			}
			else{i++;}
		}
		else{i++;}
	}
	return check;
}  //}}}




function selectionCoordinatesExpanded(){
	getSelectionCoordinates(xpoly,ypoly);
	xpoly = Array.concat(xpoly,xpoly[0]);
	ypoly = Array.concat(ypoly,ypoly[0]);
	count = xpoly.length;
	for(i=0; i<xpoly.length-1; i++){
		dist = maxOf(abs(xpoly[i]-xpoly[i+1]),abs(ypoly[i]-ypoly[i+1]));
		count += dist-1;
	}
	xpoints = newArray(count);
	ypoints = newArray(count);
	xpoints[0] = xpoly[0];
	ypoints[0] = ypoly[0];
	
	n = 1;
	for(i=1; i<xpoly.length; i++){
		dist = maxOf(abs(xpoly[i]-xpoly[i-1]),abs(ypoly[i]-ypoly[i-1]));
		if(dist>1){
			x_vector = xpoly[i]-xpoly[i-1]; // one will be 0; the other >1
			y_vector = ypoly[i]-ypoly[i-1];
			points = dist-1;
			for(j=1; j<dist; j++){
				xpoints[n] = xpoly[i-1]+j*x_vector/dist;
				ypoints[n] = ypoly[i-1]+j*y_vector/dist;
				n++;
			}
		}
		xpoints[n] = xpoly[i]; ypoints[n] = ypoly[i]; n++;
	}
	xpoints = Array.trim(xpoints,xpoints.length-1); // drop the last point
	ypoints = Array.trim(ypoints,ypoints.length-1); // drop the last point
	run("Select None");
	makeSelection("polygon", xpoints, ypoints);
	return xpoints.length;
}




/**

// Note: outside of image, pixel value is 0
function findFirst(x,y,minvalue,findvalue) {
	count = 1;
	x_list = newArray(x);
	y_list = newArray(y);
	i = 0;
	while(i<x_list.length){
		x = x_list[i]; y = y_list[i];
		if(getPixel(x,y+1)>=minvalue){count++;}
		if(getPixel(x+1,y+1)>=minvalue){count++;}
		if(getPixel(x+1,y)>=minvalue){count++;}
		if(getPixel(x+1,y-1)>=minvalue){count++;}
		if(getPixel(x,y-1)>=minvalue){count++;}
		if(getPixel(x-1,y-1)>=minvalue){count++;}
		if(getPixel(x-1,y)>=minvalue){count++;}
		if(getPixel(x-1,y+1)>=minvalue){count++;}
		i++;
	}
		
	return count;
}

**/

//}}}

// 026b LWR Functions //{{{

// Edge Interpolate
// Moves each point to half-way along each segment
// e.g. (0,1) to (0,2), the first point becomes (0,1.5)
var xpoints = newArray(0);
var ypoints = newArray(0);
function polyEdgeSelectionInterpolate(xpoints,ypoints){
	xtemp = newArray(xpoints.length);
	ytemp = newArray(xpoints.length);
	for(i=0; i<xpoints.length-1; i++){
		xtemp[i] = 0.5 * (xpoints[i] + xpoints[i+1]);
		ytemp[i] = 0.5 * (ypoints[i] + ypoints[i+1]);
	}
	xtemp[xtemp.length-1] = 0.5 * (xpoints[xtemp.length-1] + xpoints[0]);
	ytemp[ytemp.length-1] = 0.5 * (ypoints[ytemp.length-1] + ypoints[0]);
	xpoints = xtemp; ypoints = ytemp;
}

//getSelectionCoordinates(xpoints,ypoints);
//polyEdgeSelectionInterpolate(xpoints,ypoints);
//for(i=0; i<xpoints.length-1; i++){
//print(i + " " + xpoints[i] + " " + ypoints[i]);
//}
  
// Actual shortest distance between a finite line segment and a point
// (x1,y1) & (x2,y2) are the line segment (x3,y3) is the point.
// based on quano's answer: http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
function PointToLineSegment(x1,y1,x2,y2,x3,y3){
	px = x2 - x1; py = y2 - y1;
	something = px * px + py * py;
	u = ((x3 - x1) * px + (y3 - y1) * py) / something;
	if(u>1){u = 1;} else if(u<0){u = 0;}
	x = x1 + u * px; y = y1 + u * py;
	dx = x - x3; dy = y - y3;
	dist = sqrt(dx*dx + dy*dy);
	result = newArray(dist,x,y);
	return result;
}

// Slope Function
function slopeArray(x,y,nmax,maxslope){
	slope = newArray(x.length);
	for(k=1; k<x.length-1; k++){
		n = minOf(minOf(k,nmax),x.length-1-k);
		sum_xy = 0; sum_x = 0; sum_y = 0; sum_x2 = 0; count = 0;
		for(i=-n; i<=n; i++){
			sum_y += y[i+k]; sum_x += x[i+k]; sum_xy += x[i+k]*y[i+k]; sum_x2 += x[i+k]*x[i+k]; count++;
		}
		if((sum_x2 - sum_x * sum_x / count) == 0){slope[k]=maxslope;}
		else{slope[k] = (sum_xy - sum_x * sum_y / count) / (sum_x2 - sum_x * sum_x / count);}
	}
	slope[0] = slope[1];
	slope[x.length-1] = slope[x.length-2];
	return slope;
}

/**
The goal

http://stackoverflow.com/questions/4977491/determining-if-two-line-segments-intersect?lq=1  #!# this one's the best.

http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect

http://www.loria.fr/~lazard/ARC-Visi3D/Pant-project/files/Line_Segment_Line.html

**/

function segmentIntersection(xa1,ya1,xa2,ya2,xb1,yb1,xb2,yb2){
	// Note that "B" vector is the 
	// p(t) = u0 + a * v0			// "A"" vector
	// p(s) = u1 + b * v1			// "B" vector
	x00 = xa1; y00 = ya1; 			// u0 = (x00,y00)
	x01 = xa2 - xa1; y01 = ya2 - ya1;	// v0 = (x01,y01)
	x10 = xb1; y10 = yb1;			// u1 = (x10,y10)
	x11 = xb2 - xb1; y11 = yb2 - yb1;	// v1 = (x11,y11)
	crossproduct2D = x01*y11 - y01*x11;	// if crossproduct2D == 0, the lines are parallel
	//print(crossproduct2D);
	if(crossproduct2D != 0){
		d = x11 * y01 - x01 * y11;				// determinant
		b = (1/d) * ( (x00 - x10)*y01 - (y00 - y10)*x01 ); 	// indicates position on the "B" vector
		a = (1/d) * ( (x00 - x10)*y11 - (y00 - y10)*x11 ); 	// indicates position on the "A" vector
		xint = x00 + a*x01; yint = y00 + a*y01;			// the (xint,yint) where the lines intersect
		if(a>=1 && b<=1 && b>=0){correct = 1;}else{correct = 0;}// this should preclude all values except the desired intersection
		//if(correct){
		//print("correct: " + correct + " & CP: " + crossproduct2D + " & a: " + a + " & b: " + b + " & xint: " + xint + " & yint: " + yint);
		//}
		// CP, a, b, xint, yint,
		return newArray(correct,crossproduct2D,a,b,xint,yint);
	}
	else { 
		return newArray(0,0,0,0,0,0);
		//print("correct: " + 0 + " & CP: " + 0 + " & a: " + (1/0) + " & b: " + (1/0) + " & xint: " + "---" + " & yint: " + "---");
	}
	// Note: if I was solely interested in overlapping segments, could use box/range overlap to test...
}

// g = global... the original array
// rerun = true if it needs to be rerun
function MinimaNotFoundRerunNecessary(sub_start,sub_end,g_entry,g_start,g_end){
	if(g_entry==0 || g_entry==g_end-1){return false;} // no need to rerun; probably an end point
	else if(g_entry-sub_start>0 && g_entry-sub_start<(sub_end-sub_start)-1){return false;} // no need to rerun; it's in the middle somewhere
	else{return true;} // need to rerun
}

// Cross Product
// When u = (x_edge-x_skel, y_edge-y_skel) & v = (1, slope)
function nCrossProduv(ux,uy,vx,vy){
	cross_product = ux*vy - vx*uy;
	if(cross_product==0){sign_CP = 0;}else{sign_CP = cross_product/abs(cross_product);}
	return sign_CP;
}


//}}}

// 027 Downsample Functions //{{{

// Fibonacci Downsample:
function downsampleFibonacci(array,factor){
	new_length = floor(array.length / factor);
	temp = newArray(new_length);
	density = 1 / factor; df = 1;
	fibonacci = newArray(1,1); f_sum = 2;
	while(df>=density){
		fibonacci = Array.concat(fibonacci,fibonacci[fibonacci.length-1]+fibonacci[fibonacci.length-2]);
		f_sum += fibonacci[fibonacci.length-1];
		df = fibonacci.length / f_sum;
		//df = fibonacci.length / fibonacci[fibonacci.length-1];
	}
	f_sum = f_sum - fibonacci[fibonacci.length-1];
	fibonacci[fibonacci.length-1] = floor(factor*fibonacci.length) - f_sum - 1;
	f_sum += fibonacci[fibonacci.length-1];
	//df = fibonacci.length / f_sum;
	//Array.print(fibonacci);
	fi = 0; fii = 0;
	for(i=0; i<new_length; i++){
		fii = fibonacci[fi] + fii; if(fii>array.length-1){fii=array.length-1;}
		temp[i] = array[fii]; //print(i + " " + temp[i]);
		if(fi<fibonacci.length-1){fi++;}else{fi=0;}
	}
	return temp;
}

// Random Selector (for random downsampling)
function random_selector(factor,initial_points){
	ds_points = round(initial_points / factor * 2);
	downsampled_points = newArray(ds_points);
	sum = 0; i = 0;
	while(sum<initial_points){
		item = round(2*factor*random);
		while(item==0){ item = round(2*factor*random); }
		sum += item;
		downsampled_points[i] = item;
		i++;
	}
	downsampled_points = Array.trim(downsampled_points,i-1);
	ds_result = newArray(downsampled_points.length);
	ds_result[0] = downsampled_points[0];
	for(i=1; i<downsampled_points.length; i++){
		ds_result[i] = downsampled_points[i] + ds_result[i-1];
	}
	return ds_result;
}

// Linear Downsample (every Nth item is taken)
/** downsample function is imperfect... misses last value **/
function downsample(array,factor){
	new_length = floor(array.length / factor);
	print(new_length);
	temp = newArray(new_length);
	for(i=0; i<new_length; i++){
		temp[i] = array[i*factor];
	}
	return temp;
}

// 027 //}}}

// 028 Base-62 Converter //{{{
// Base 62 Converter (new)

function base62converter(number,convert){
	c = newArray("0","1","2","3","4","5","6","7","8","9","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z");
	base = c.length;
	if(convert){
		on = true; basenumber = "";
		while(on){
			num = floor(number/base);
			rem = number - num*base;
			basenumber = c[rem] + basenumber;
			if(num > 0){ number = num;} else {on = false;}
		}
	}
	else{
		on = true; ibase = 0; basenumber = 0;
		while(on){
			digit = substring(number,lengthOf(number)-1);
			number = substring(number,0,lengthOf(number)-1);
			for(i=0; i<base; i++){
				if(digit==c[i]){
					if(matches(digit,"[0-9]") && matches(c[i],"[0-9]")){ basenumber += i*pow(base,ibase); ibase++;}
					if(matches(digit,"[a-z]") && matches(c[i],"[a-z]")){ basenumber += i*pow(base,ibase); ibase++;}
					if(matches(digit,"[A-Z]") && matches(c[i],"[A-Z]")){ basenumber += i*pow(base,ibase); ibase++;}
				}
			}
			if(lengthOf(number)==0){on = false;}
		}
	}
	return basenumber;
}
// 028 //}}}

// 029 ArraySmooth (Median & Window Filters) //{{{

function arraySmoothWindowSymmetric(histogram,iterations,plusminus) {
	new_histogram = Array.copy(histogram);
	old_histogram = Array.copy(histogram);
	count = 0;
	while(count<iterations){
		for(n=0; n<histogram.length; n++){
			left_limit = maxOf(0,n-plusminus); right_limit = minOf(histogram.length-1,n+plusminus);
			pm = minOf(abs(left_limit-n),abs(right_limit-n));
			left_limit = maxOf(0,n-pm); right_limit = minOf(histogram.length-1,n+pm);
			sum = 0;
			for(m=left_limit; m<right_limit+1; m++){ sum += old_histogram[m]; }
			new_histogram[n] = sum/(right_limit-left_limit+1);
		}
		
		old_histogram = Array.copy(new_histogram);
		count++;
	}
	return new_histogram;
}


function arraySmoothFilterSymmetric(histogram,iterations,plusminus) {
	new_histogram = Array.copy(histogram);
	old_histogram = Array.copy(histogram);
	count = 0;
	while(count<iterations){
		for(n=0; n<histogram.length; n++){
			left_limit = maxOf(0,n-plusminus); right_limit = minOf(histogram.length-1,n+plusminus);
			pm = minOf(abs(left_limit-n),abs(right_limit-n));
			left_limit = maxOf(0,n-pm); right_limit = minOf(histogram.length-1,n+pm);
			items = newArray(right_limit-left_limit+1); //sum = 0;
			for(m=left_limit; m<right_limit+1; m++){ items[m-left_limit] = old_histogram[m]; }
			items_sorted = Array.sort(items);
			median_item = items_sorted[pm];
			new_histogram[n] = median_item; //sum/(right_limit-left_limit+1);
		}
		
		old_histogram = Array.copy(new_histogram);
		count++;
	}
	return new_histogram;
}


function arraySmoothWindow(histogram,iterations,plusminus) {
	new_histogram = Array.copy(histogram);
	old_histogram = Array.copy(histogram);
	count = 0;
	while(count<iterations){
		for(n=0; n<histogram.length; n++){
			left_limit = maxOf(0,n-plusminus); right_limit = minOf(histogram.length-1,n+plusminus);
			sum = 0;
			for(m=left_limit; m<right_limit+1; m++){ sum += old_histogram[m]; }
			new_histogram[n] = sum/(right_limit-left_limit+1);
		}
		
		old_histogram = Array.copy(new_histogram);
		count++;
	}
	return new_histogram;
}

function arraySmoothFilter(histogram,iterations,plusminus) {
	new_histogram = Array.copy(histogram);
	old_histogram = Array.copy(histogram);
	count = 0;
	while(count<iterations){
		for(n=0; n<histogram.length; n++){
			left_limit = maxOf(0,n-plusminus); right_limit = minOf(histogram.length-1,n+plusminus);
			items = newArray(right_limit-left_limit+1); //sum = 0;
			for(m=left_limit; m<right_limit+1; m++){ items[m-left_limit] = old_histogram[m]; }
			items_sorted = Array.sort(items);
			median_item = items_sorted[plusminus];
			new_histogram[n] = median_item; //sum/(right_limit-left_limit+1);
		}
		
		old_histogram = Array.copy(new_histogram);
		count++;
	}
	return new_histogram;
}

// 029 //}}}

// 030 Object Histogram Stuff //{{{

var xObject = newArray(0);
var yObject = newArray(0);

var eightx = newArray(0,1,1,1,0,-1,-1,-1);
var eighty = newArray(1,1,0,-1,-1,-1,0,1);

function objectPixels(xstart,ystart,tempvalue){
	/** Try using index instead of concat **/
	w = getWidth(); h = getHeight();
	xObject = newArray(w*h); // clear to ensure empty
	yObject = newArray(w*h); // clear to ensure empty
	xObject[0] = xstart; 
	yObject[0] = ystart;
	value = getPixel(xstart,ystart);
	setForegroundIndex(tempvalue);
	floodFill(xstart,ystart,"8-connected");
	setPixel(xstart,ystart,value);
	//xtemp = newArray(8);
	//ytemp = newArray(8);
	inOperation = true; index = 0; appendindex = 1;
	while(inOperation){
		x = xObject[index]; y = yObject[index];
		//count = 0;
		for(i=0; i<8; i++){
			if(getPixel(x+eightx[i],y+eighty[i])==tempvalue){
				xObject[appendindex] = x+eightx[i];
				yObject[appendindex] = y+eighty[i];
				appendindex++;
				setPixel(x+eightx[i],y+eighty[i],value);
			}
		}
		if(index<appendindex-1){index++;}else{inOperation = false;}
	}
	print(xObject.length);
	xObject = Array.trim(xObject,appendindex);
	yObject = Array.trim(yObject,appendindex);
	return xObject.length;
}

function objectHistogram(xstart,ystart,tempvalue,binary_image,original_image){
	original = getImageID();
	selectImage(binary_image);
	count = objectPixels(xstart,ystart,tempvalue);
	values_array = newArray(count);
	histogram_array = newArray(256);
	selectImage(original_image);
	for(i=0; i<count; i++){
		values_array[i] = getPixel(xObject[i],yObject[i]);
	}
	for(i=0; i<count; i++){
		histogram_array[values_array[i]] += 1;
	}
	selectImage(original);
	return histogram_array;
}

function histogramRangeCount(first,last,histogram){
	count = 0;
	sum = 0;
	for(i=0;i<256;i++){
		if(i>=first && i<=last){ count += histogram[i]; }
		sum += histogram[i];
	}
	proportion = count / sum;
	results = newArray(proportion,count,sum);
	return results;
}

//}}}

// 031 OPA Specific: Components & followTwoEncodeVS (variant) //{{{

function components(array,limit){
	if(array.length < limit){
		values = newArray(array.length);
		value = array[array.length-1] - array[0];
		Array.fill(values,value);
	}
	else{
		start = round(limit/2);
		end = limit-start-1;
		values = newArray(array.length);
		for(i=start; i<array.length-end; i++){
			values[i] = array[i+end] - array[i-start];
		}
		for(i=0; i<start; i++){
			values[i] = values[start];
		}
		for(i=array.length-end; i<array.length; i++){
			values[i] = values[array.length-end-1];
		}
	}
	return values;
}

// this function has been modified
function followTwoEncodeVS(x,y,value,valueswap) {
	condition = true; pixels = 0; u = x; v = y; a = 0; b = 0;
	encoded = newArray(0);
	while(condition){
		// Find which direction to go:
		count = 0; i = 0; j = 0; ip = u-x ; jp = v-y ; 
		// count = neighbours; i,j = direction to be moved; u,v = previous position; a = t-movement; b = x-movemnt;
		if(getPixel(x,y+1)>=value){count++; if(ip == 0 && jp == 1){ }else{ i = 0; j = 1; }}
		if(getPixel(x+1,y+1)>=value){count++; if(ip == 1 && jp == 1){ }else{ i = 1; j = 1; }}
		if(getPixel(x+1,y)>=value){count++; if(ip == 1 && jp == 0){ }else{ i = 1; j = 0; }}
		if(getPixel(x+1,y-1)>=value){count++; if(ip == 1 && jp == -1){ }else{ i = 1; j = -1; }}
		if(getPixel(x,y-1)>=value){count++; if(ip == 0 && jp == -1){ }else{ i = 0; j = -1; }}
		if(getPixel(x-1,y-1)>=value){count++; if(ip == -1 && jp == -1){ }else{ i = -1; j = -1; }}
		if(getPixel(x-1,y)>=value){count++; if(ip == -1 && jp == 0){ }else{ i = -1; j = 0; }}
		if(getPixel(x-1,y+1)>=value){count++; if(ip == -1 && jp == 1){ }else{ i = -1; j = 1; }}
		// setPixel(x,y,value); // print(steps + " " + x + " " + y);
		// Check to see if we should continue:
		encoded = Array.concat(encoded,xy_coder(x,y,0,"enc")); if(valueswap!=value){setPixel(x,y,valueswap);}
		if(pixels > 0 && count == 2){ u = x; v = y; x += i; y+= j; pixels++;} 
		else if(pixels==0 && count<=2){ u = x; v = y; x += i; y+= j; pixels++;}
		else {condition = false;  pixels++; i = 0; j = 0; }
		if(abs(i)+abs(j) == 1){a++;} else if(abs(i)+abs(j) == 2){b++;}

	}
	distance = a + sqrt(2)*b; result = newArray(x,y,distance,pixels);
	return encoded;
}

// 030 //}}}

// 032 AF Particle Analysis Excerpt additional functions {{{

/** New Functions **/

function furthestNonEdgePixel(x,y,w,h,cvalue){
	//if(w-1-x < x){xedge = w;} else {xedge = 0;}
	//if(h-1-y < y){yedge = h;} else {yedge = 0;}
	doWand(x,y); getSelectionCoordinates(xvalues,yvalues);
	run("Select None");
	xv = x; yv = y;
	for(sc=0; sc<xvalues.length; sc++){
		//dist = minOf(abs(xvalues[sc]-xedge),abs(yvalues[sc]-yedge));
		dist = minOf(minOf(abs(xvalues[sc]-w-1),abs(yvalues[sc]-h-1)),minOf(xvalues[sc],yvalues[sc]));
		if(sc==0){mdist = dist;}
		if(dist>mdist){xv = xvalues[sc]; yv = yvalues[sc]; mdist = dist;}
	}
	if(getPixel(xv,yv)==cvalue){}
	else if(getPixel(xv+1,yv)==cvalue){xv = xv+1;}
	else if(getPixel(xv-1,yv)==cvalue){xv = xv-1;}
	else if(getPixel(xv,yv+1)==cvalue){yv = yv+1;}
	else if(getPixel(xv,yv-1)==cvalue){yv = yv-1;}
	else if(getPixel(xv+1,yv+1)==cvalue){xv = xv+1; yv = yv+1;}
	else if(getPixel(xv-1,yv-1)==cvalue){xv = xv-1; yv = yv-1;}
	else if(getPixel(xv-1,yv+1)==cvalue){xv = xv-1; yv = yv+1;}
	else if(getPixel(xv+1,yv-1)==cvalue){xv = xv+1; yv = yv-1;}
	return newArray(xv,yv);
}

function arrayFindValue(array,value){
	i = 0;
	while(i<array.length){
		if(array[i]==value){ivalue = i; return ivalue;}
		i++;
	}
	return -1;
}


function imagesaver(switch,image_ID,directory,filename,filetype){
	if(switch){
		if(filetype = "tif"){ ext = "tif";} else { ext = filetype; }
		selectImage(image_ID); temp_title = getTitle();	saveAs(filetype,directory+filename+"." + ext); rename(temp_title);
	}
}

//}}}

// 033 Defect Drawing Tools & Shapes //{{{

// SHAPES
// All shapes are normalized, such that branches are ~ 1 px wide; they should be scaled according to the linewidth
// See: Special_Defect_Markers_03.ijm for development / tests

// TREFOIL for 3-Junctions
Trefoil_x = newArray(0.0006, -0.0921, -0.1786, -0.2569, -0.3252, -0.3815, -0.4240, -0.4509, -0.4603, -0.4603, -0.9212, -0.9969, -1.0583, -1.1049, -1.1358, -1.1505, -1.1482, -1.1283, -1.0900, -1.0355, -0.9689, -0.8929, -0.8100, -0.7228, -0.6337, -0.5453, -0.4603, 0.0006, 0.4600, 0.5451, 0.6334, 0.7225, 0.8098, 0.8927, 0.9687, 1.0352, 1.0897, 1.1280, 1.1479, 1.1502, 1.1356, 1.1046, 1.0581, 0.9966, 0.9210, 0.4600, 0.4600, 0.4507, 0.4240, 0.3817, 0.3257, 0.2577, 0.1797, 0.0934);
Trefoil_y = newArray(-1.2579, -1.2486, -1.2217, -1.1791, -1.1228, -1.0546, -0.9763, -0.8898, -0.7970, -0.2657, 0.0014, 0.0558, 0.1220, 0.1975, 0.2798, 0.3665, 0.4551, 0.5431, 0.6280, 0.7037, 0.7651, 0.8117, 0.8426, 0.8573, 0.8550, 0.8350, 0.7967, 0.5311, 0.7967, 0.8351, 0.8552, 0.8578, 0.8434, 0.8127, 0.7664, 0.7052, 0.6296, 0.5445, 0.4562, 0.3671, 0.2798, 0.1969, 0.1209, 0.0544, -0.0001, -0.2657, -0.7970, -0.8898, -0.9763, -1.0546, -1.1228, -1.1791, -1.2217, -1.2486);

// PLUS for 4-junctions
Plus_x = newArray(0, -0.1008, -0.1946, -0.2796, -0.3536, -0.4146, -0.4607, -0.4898, -0.5000, -0.5000, -1.1000, -1.2008, -1.2946, -1.3796, -1.4536, -1.5146, -1.5607, -1.5898, -1.6000, -1.5898, -1.5607, -1.5146, -1.4536, -1.3796, -1.2946, -1.2008, -1.1000, -0.5000, -0.5000, -0.4898, -0.4607, -0.4146, -0.3536, -0.2796, -0.1946, -0.1008, 1.000E-6, 0.1008, 0.1946, 0.2796, 0.3536, 0.4146, 0.4607, 0.4898, 0.5000, 0.5000, 1.1000, 1.2008, 1.2946, 1.3796, 1.4536, 1.5146, 1.5607, 1.5898, 1.6000, 1.5898, 1.5607, 1.5146, 1.4536, 1.3796, 1.2946, 1.2008, 1.1000, 0.5000, 0.5000, 0.4898, 0.4607, 0.4146, 0.3536, 0.2796, 0.1946, 0.1008);
Plus_y = newArray(-1.6000, -1.5900, -1.5610, -1.5150, -1.4540, -1.3800, -1.2950, -1.2010, -1.1002, -0.5002, -0.5002, -0.4902, -0.4612, -0.4152, -0.3542, -0.2802, -0.1952, -0.1012, -0.0004, 0.1004, 0.1944, 0.2794, 0.3534, 0.4144, 0.4604, 0.4894, 0.4994, 0.4994, 1.0994, 1.2002, 1.2942, 1.3792, 1.4532, 1.5142, 1.5602, 1.5892, 1.5992, 1.5892, 1.5602, 1.5142, 1.4532, 1.3792, 1.2942, 1.2002, 1.0994, 0.4994, 0.4994, 0.4894, 0.4604, 0.4144, 0.3534, 0.2794, 0.1944, 0.1004, -0.0004, -0.1012, -0.1952, -0.2802, -0.3542, -0.4152, -0.4612, -0.4902, -0.5002, -0.5002, -1.1002, -1.2010, -1.2950, -1.3800, -1.4540, -1.5150, -1.5610, -1.5900);

// STAR for 5-junctions and above
Star_x = newArray(0.7478, 1.2993, 1.9970, 2.6314, 2.8596, 2.9934, 2.9848, 2.8308, 2.2545, 1.6000, 1.2031, 1.2934, 1.6543, 1.9521, 1.9730, 1.8529, 1.6448, 1.4296, 0.9795, 0.5052, 0.0096, -0.4862, -0.9608, -1.4112, -1.6264, -1.8345, -1.9545, -1.9334, -1.6353, -1.2740, -1.1834, -1.5801, -2.2343, -2.8105, -2.9644, -2.9729, -2.8390, -2.6108, -1.9764, -1.2788, -0.7272, -0.4766, -0.4064, -0.3120, -0.1919, 0.0109, 0.2137, 0.3337, 0.4276, 0.4974);
Star_y = newArray(-1.0059, -0.9216, -1.0744, -1.1814, -1.1292, -0.9595, -0.7434, -0.5672, -0.2812, 0.0048, 0.3970, 0.9476, 1.5639, 2.1343, 2.3674, 2.5471, 2.6057, 2.5136, 2.0539, 1.5199, 1.2636, 1.5196, 2.0533, 2.5128, 2.6047, 2.5460, 2.3662, 2.1331, 1.5629, 0.9469, 0.3963, 0.0039, -0.2825, -0.5689, -0.7452, -0.9613, -1.1309, -1.1830, -1.0756, -0.9224, -1.0063, -1.5049, -2.2156, -2.8521, -3.0529, -3.1278, -3.0528, -2.8519, -2.2153, -1.5046);


// Function: Rotate, translate, scale, & select combation
// note: selection must already be centred
// half_true is the condition for the 0.5 offset.
function selectionCombo(xPoints,yPoints,angle,degrees_true,scale_factor,x_offset,y_offset,half_true){
	x_centre = 0; y_centre = 0;
	if(half_true){half = 0.5;}else{half = 0;}
	if(degrees_true){ angle = PI*angle/180; }
	cos_theta = cos(angle);
	sin_theta = sin(angle);
	
	x_new = Array.copy(xPoints);
	y_new = Array.copy(yPoints);
	
	for(i=0; i<x_new.length; i++){
		x_new[i] = scale_factor*(cos_theta*(xPoints[i]-x_centre) - sin_theta*(yPoints[i]-y_centre)) + x_centre + x_offset + half;
		y_new[i] = scale_factor*(sin_theta*(xPoints[i]-x_centre) + cos_theta*(yPoints[i]-y_centre)) + x_centre + y_offset + half;
	}
	makeSelection("polygon",x_new,y_new);
}

//}}}

//// END OF FUNCTIONS //}}}



// USER OPTIONS  //{{{

	// Question: Batch Process //{{{
		options_single_multi = newArray("Single","Multi");
		Dialog.create("Single or Multi Mode");
			Dialog.addMessage("Multi Mode will analyze a folder of TIF or PNG images. \n Single mode will only process one image. \n The image must already be open & active.");
			Dialog.addChoice("Choose:",options_single_multi);
			//Dialog.addMessage("Batch Mode. (Hide Images)");
			Dialog.addCheckbox("Batch Mode", true);
			Dialog.show();
		single_or_multi = Dialog.getChoice();
		batchmode_choice = Dialog.getCheckbox(); //}}}
		//!@#$%

	// Question: Resolution //{{{
  		options_resolution = newArray("NINT_s4800","Embedded","Specified");
  		options_resolution_units = newArray("nm",getInfo("micrometer.abbreviation"));
  		default_value_resolution = 0.500; // nm per pixel;
  		Dialog.create("Determine Image Resolution");
			Dialog.addMessage("How will the image resolution be determined? \nIf 'Specified' is selected, enter a value and \nchoose a unit.");
			Dialog.addChoice("Choose Method:",options_resolution);
			Dialog.addNumber("1 pixel equals:",default_value_resolution,4,9,"units");
			Dialog.addChoice("Unit:",options_resolution_units);
			Dialog.show();
		resolution_method = Dialog.getChoice();
		resolution_specified = Dialog.getNumber();
		resolution_units = Dialog.getChoice(); //}}}
		//!@#$%
		
	// Question: Wire Period & Method //{{{
		options_period = newArray("Automatic","WLS","FFT","LTL","Specified");
		options_period_units = newArray("pixels","nm",getInfo("micrometer.abbreviation"));
		Dialog.create("Determine Line Period & Width");
			Dialog.addMessage("Automatic = WLS + FFT \nWSL = Weighted Least Squares \nFFT = Fast Fourier Transform \nLTL = Line-to-Line \nSpecified = user-input period");
			Dialog.addChoice("Method:",options_period);
			Dialog.addNumber("Period:",10,4,12,"units");
			Dialog.addChoice("Unit:",options_period_units);
			Dialog.addNumber("Range: ±",0,4,6,"units");
			Dialog.addMessage("Note: if Range > 0, the specified \nperiod will be utilized for 'Automatic'");
			Dialog.addMessage("Defined variables allow for periods \nranging from "+period_range_min_nm+" nm to "+period_range_max_nm+" nm.\nRequires code edit to change.");
			Dialog.show();
		period_method = Dialog.getChoice();
		period_specified = Dialog.getNumber();
		period_units = Dialog.getChoice();
		period_range = Dialog.getNumber();  //}}}
		//!@#$%
		
	// Question: Smoothing & Noise Reduction //{{{
		options_smoothing = newArray("Automatic","None","Median","Gaussian");
		Dialog.create("Additional Smoothing Procedures");
			Dialog.addMessage("Choose what smoothing, if any,\nis required for the image.\nFor Median or Gaussian,\n choose a radius.");
			Dialog.addChoice("Method:",options_smoothing);
			Dialog.addNumber("Radius:",5,1,4,"pixels");
			Dialog.addCheckbox("Relative", true);
			Dialog.addNumber("Radius:",20,1,4,"% of period");
			Dialog.show();
		smoothing_method = Dialog.getChoice();
		smoothing_radius = Dialog.getNumber();
		smoothing_relative_option = Dialog.getCheckbox();
		smoothing_relative_percent = Dialog.getNumber(); //}}}
		
	// Question: Thresholding //{{{
		options_thresholding = newArray("Auto-Local","Custom","Built-in");
		ij_thresholding_options = newArray("Default","Huang","Intermodes","IsoData","IJ_IsoData","Li","Mean","Minimum","Otsu","Triangle");
		auto_local_options = newArray("Otsu","Bernsen","Contrast","Mean","Median","MidGrey","Niblack","Phansalkar","Sauvola");
		auto_local_p1 = 0; auto_local_p2 = 0; auto_local_radius = 15;
		Dialog.create("Choose Thresholding Method");
			Dialog.addMessage("Choose the thresholding method to apply,\nto create the binary image.\n");
			Dialog.addChoice("Method Type:",options_thresholding);
			Dialog.addChoice("Method:",ij_thresholding_options);
			Dialog.addMessage("Auto-Local options:");
			Dialog.addChoice("Auto-Local Type:",auto_local_options);
			Dialog.addCheckbox("Auto parameters", true);
			Dialog.addNumber("Radius:",auto_local_radius,2,4,"pixels");
			Dialog.addNumber("Parameter 1:",auto_local_p1,3,4,"");
			Dialog.addNumber("Parameter 2:",auto_local_p2,3,4,"");
			Dialog.addMessage("Override the above options &\n use radius=1.5*period.\n");
			Dialog.addCheckbox("Use relative instead", true); //option simply uses r = 1.5 * period.
			Dialog.show();
		thresholding_choice = Dialog.getChoice();
		ij_thresholding_choice = Dialog.getChoice(); 
		auto_local_choice = Dialog.getChoice();
		auto_local_auto = Dialog.getCheckbox();
		auto_local_radius = Dialog.getNumber();
		auto_local_p1 = Dialog.getNumber();
		auto_local_p2 = Dialog.getNumber();
		auto_local_relative = Dialog.getCheckbox();//}}}
		
	// Question: Aditional Analyses //{{{
		Dialog.create("Choose Additional Analyses");
			Dialog.addMessage("Choose whether the program will \nperform any addtional analyses.\nNote that this will increase run time.");
			Dialog.addCheckbox("Orientational Domain Map", true);
			Dialog.show();
		domain_mapping = Dialog.getCheckbox();
	//}}}
	
	// Question: Cropping //{{{
		crop_message = "";
		while(crop_message != " "){
		cropping_options = newArray("No","Yes");
		cropping_step_opt = newArray("Stage 1","Stage 2");
		cropping_unit_opt = newArray("Percent","Pixels");
		cropping_uniform_opt = newArray("Uniform","Non-Uniform");
		Dialog.create("Crop Options");
			Dialog.addMessage(crop_message + "Automatic cropping of images.\nNote: for NINT images,\ncrop is automatic.\n.");
			Dialog.addChoice("Crop:",cropping_options);
			Dialog.addChoice("When:",cropping_step_opt); 
			Dialog.addMessage("Percent interpreted as % value.\nPixels interpreted as px value.\n.");
			Dialog.addChoice("% or px:",cropping_unit_opt);
			Dialog.addChoice("Uniform:",cropping_uniform_opt);
			Dialog.addNumber("Crop each edge:",10,2,10,"(% or pixels)");
			Dialog.show();
		crop_choice = Dialog.getChoice();
		crop_step = Dialog.getChoice();
		crop_unit = Dialog.getChoice();
		crop_uniform = Dialog.getChoice();
		crop_value = Dialog.getNumber();
		if(crop_choice == "Yes"){
			if(crop_uniform=="Non-Uniform"){
				Dialog.create("How to crop each edge?");
					Dialog.addNumber("X: Left edge",crop_value,2,10,crop_unit);
					Dialog.addNumber("X: Right edge",crop_value,2,10,crop_unit);
					Dialog.addNumber("Y: Top edge",crop_value,2,10,crop_unit);
					Dialog.addNumber("Y: Bottom edge",crop_value,2,10,crop_unit);
					Dialog.show();
				crop_x_left = Dialog.getNumber();
				crop_x_right = Dialog.getNumber();
				crop_y_top = Dialog.getNumber();
				crop_y_bottom = Dialog.getNumber();
				if(crop_unit == "Percent"){ if((crop_x_left + crop_x_right >= 50) || (crop_y_top + crop_y_bottom >= 50)){ crop_message = "Error. Retry.\n"; } else{ crop_message = " "; }} else{ crop_message = " "; }
			}	
			else if(crop_unit == "Percent"){ if(crop_value >= 50){ crop_message = "Error. Retry.\n"; } else{ crop_message = " "; }} else{ crop_message = " "; }
		}
		else { crop_message = " "; }
		}
		//}}}

		
  
  // Image Lists  //{{{
  	if(single_or_multi=="Single"){
  		if(nImages() >= 1){
				image_choice = getImageID();
				image_list = newArray(1);
				image_list[0] = getTitle();
  		} else { exit ("No images are currently open! \n Open an image to be analyzed, \n then restart the macro."); }
  	}
		else if(single_or_multi=="Multi"){
			// List All images in the directory... 
			image_directory = getDirectory("Choose a Folder");
			file_list = getFileList(image_directory);
			image_count = 0;
			for(i=0; i<file_list.length; i++){
				if(endsWith(toLowerCase(file_list[i]),"tif")||endsWith(toLowerCase(file_list[i]),"png")){
					image_count++;
				}
			}
			image_list = newArray(image_count);
			image_count = 0;
			for(i=0; i<file_list.length; i++){
				if(endsWith(toLowerCase(file_list[i]),"tif")||endsWith(toLowerCase(file_list[i]),"png")){
					image_list[image_count] = file_list[i];
					image_count++;
				}
			}
		}
		
		
  // END OF Image Lists //}}}
  
  // Output Folder //{{{
  	// Folder will be named according to: Output_YYYYMMDD_HHMMSS
		getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
			year_dd = toString(year); month++;
			if(month<10){month_dd ="0" + toString(month);} else {month_dd = toString(month);}
			if(dayOfMonth<10){day_dd ="0" + toString(dayOfMonth);} else {day_dd = toString(dayOfMonth);}
			if(hour<10){hour_dd ="0" + toString(hour);} else {hour_dd = toString(hour);}
			if(minute<10){minute_dd ="0" + toString(minute);} else {minute_dd = toString(minute);}
			if(second<10){second_dd ="0" + toString(second);} else {second_dd = toString(second);}
			time_stamp = year_dd + month_dd + day_dd + "_" + hour_dd + minute_dd + second_dd;
		output_folder = "Output_" + time_stamp;
		save_location = getDirectory("Choose save location for Output folder.");
		File.makeDirectory(save_location + output_folder);
		save_folder = save_location + output_folder + File.separator;  //}}}
		
		// LOG FILE  //{{{
		print("\\Clear");  // Clear the Log Window
		print(program_name + " " + program_version + " " + modification_date);
		// Time-stamp: gets everything exactly as I want it.
		print("Time-stamp: "+year+"."+month_dd+"."+day_dd+" "+hour+":"+minute_dd);
		print(" ");
		// Record User Choices  //!@#$%
		

		// Print a list of images to be analyzed
		for(i=0; i<image_list.length; i++){
			print(i+": "+image_list[i]);
		}
				
		// Save Log File
		selectWindow("Log");
		saveAs("Text", save_folder + output_folder + ".txt");
		print("\\Clear");  // Clear the Log Window
		// END OF LOG FILE  //}}}
		

//// END OF USER OPTIONS //}}}

if(batchmode_choice){setBatchMode(true);}

// Multi-Image Data Storage //{{{
if(image_list.length>1){
	/** Insert Arrays Here */
	MULTI_defect_pair_density_um = newArray(image_list.length);
	MULTI_time = newArray(image_list.length);
}
// End of M.I.D.S //}}}

// ###### MULTI-IMAGE LOOP START ######  //{{{
time_zero = getTime();
for(img_i=0; img_i<image_list.length; img_i++){
	time_image_start = getTime() - time_zero; run("Clear Results"); //!@#$ not sure if absolutely necessary here, but prob a good idea
	// IMAGE IMPORT  //{{{
	if(single_or_multi == "Multi"){
		open(image_directory + image_list[img_i]);
		image000 = getImageID();
	}
	else if(single_or_multi == "Single"){
		image000 = image_choice;
	}
	series_number = img_i;
	// END OF IMAGE IMPORT  //}}}
	
	// IMAGE VIABILITY  //{{{
		viable_image = true; 
		if(bitDepth() != 8){viable_image = false;} else { run("8-bit"); } // Fix to turn "8-bit Color" images to plain "8-bit"
		//else if(true){viable_image = false;}  // Other reasons...
	// END OF IMAGE VIABILITY  //}}}
	
	// ALGORITHM CODE: VIABLE IMAGES //{{{
	if(viable_image){
			
		// DATA STRUCTURES //{{{
			output_data = newArray(0);
			output_tags = newArray(0);
			output_labels = newArray(0);
			//OUTPUT DATA + TAGS
			outputTD("Image_Number",img_i); outputTD("Image_Title.String",image_list[img_i]); outputTD("Output_Folder",output_folder);
		// DATA STRUCTURES //}}}
	
		// CREATE SUBFOLDER (for saving files) //{{{
		File.makeDirectory(save_folder + series_number);
		save_subfolder = save_folder + series_number + File.separator;  
		selectImage(image000); saveAs("tiff",save_subfolder+"image000"+".tif"); //}}}
		
	// NEW IMAGE LOG  //{{{
	print("\\Clear");
	print(series_number + ": " + image_list[series_number]);
	print("Folder: " + output_folder);
	print("Time elapsed: " + time_image_start + " msec" + "\n");
	print(program_name + " " + program_version + " " + modification_date);
	print("Time-stamp: "+year+"."+month_dd+"."+day_dd+" "+hour+":"+minute_dd + "\n");
	//OUTPUT DATA + TAGS
	output_date = year*10000 + month_dd*100 + day_dd; output_time = hour*100 + minute_dd;
	outputTD("Date",output_date); 
	outputTD("Time",output_time);
	outputTD("Version",prog_version); 
	// END OF NEW IMAGE LOG  //}}}

	// DETERMINE RESOLUTION  //{{{
	if(resolution_method == "Embedded"){
		getPixelSize(unit, pixelWidth, pixelHeight);  
		if(unit == "nm"){ nm_per_pixel = pixelWidth;}
		else if(unit == getInfo("micrometer.abbreviation") || unit == "um"){
			nm_per_pixel = 1000*pixelWidth; //!@#$ Assumes pxWidth = pxHeight
		}
		else exit("Image Resolution is not embedded in " + image_list[series_number]); 
		run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000 frame=[0 sec] origin=0,0");
	}
	else if(resolution_method == "Specified"){
		if(resolution_units == "nm"){nm_per_pixel = resolution_specified;}
		else if(resolution_units == getInfo("micrometer.abbreviation")){nm_per_pixel = 1000*resolution_specified;}
		else exit("Something's wrong here");
		run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000 frame=[0 sec] origin=0,0");
	}
	else if(resolution_method == "NINT_s4800"){
		
	// Measure Scale //{{{
	// Identifying Characters for the Scalebar
	// Modify Settings
		run("Properties...", "channels=1 slices=1 frames=1 unit=pixel pixel_width=1 pixel_height=1 voxel_depth=1.0000000 frame=[0 sec] origin=0,0");
	
	// Character Identification Settings  //{{{
	// Digits (0-9), Letters (num), & Punctuation (.) according to Particle Size Analysis
	character_array = newArray("0","1","2","3","4","5","6","7","8","9",".","n","u","m");
		
		w = getWidth(); // Determine the Image Size; different font sizes.
		if(w==1280){
			character_area_array = newArray(190,105,180,172,177,188,210,188,216,205,9,144,141,228);  //1280 version
		} else if(w==2560){
			character_area_array = newArray(876,475,810,757,814,847,947,563,941,926,49,656,654,1010);  //2560 version
		} else {
			exit("unknown image size");  // Error message
		}  //END OF Character Identification Settings //}}}
	
	// Character Recognition Routine //{{{
		if(w==1280){
			makeRectangle(1125, 912, 155, 48); // needs to change depending on image size
		} else if(w==2560){
			makeRectangle(2230, 1824, 330, 96);
		}
		// Copy ScaleBar Text & Analyse Character Sizes
		run("Duplicate...", "title=scalebar_text"); run("Invert LUT"); scalebar = getImageID();
		run("Set Measurements...", "area mean bounding redirect=None decimal=3"); 
		run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display clear record");
		selectImage(scalebar); close(); run("Select None");
		
		character_x_starts = newArray(nResults);
		character_particle_area = newArray(nResults);
		character_sequenced = newArray(nResults);
		scalebar_number_string = "";
		scalebar_unit_string = "";
		
		for(x=0; x<nResults; x++){
				character_x_starts[x] = getResult("XStart",x);
		}
		character_x_sequence = Array.rankPositions(character_x_starts); // Sort them by XStart to give Left-to-Right Sequence
		// Character Particle Sizes in Correct Sequence
		for(x=0; x<nResults; x++){
				n = character_x_sequence[x];
				character_particle_area[x] = getResult("Area",n);
		}
		// Convert Sizes to Characters(Numbers,Letters,Punctuation)
		for(n=0; n<nResults; n++){
			m = 0;
			while(m<character_area_array.length){
				if(character_particle_area[n] == character_area_array[m]){
					character_index = m;
					m = character_area_array.length;
				}
				m++;
			}
			character_sequenced[n] = character_array[character_index];
			// Separate Numbers from Units
			if(character_index<=10){
				scalebar_number_string = scalebar_number_string + character_array[character_index];
			} else {
				scalebar_unit_string = scalebar_unit_string + character_array[character_index];
			}
		}
		run("Clear Results"); selectWindow("Results"); run("Close"); 
			
		print(scalebar_number_string + " " + scalebar_unit_string);
		scalebar_number_int = parseInt(scalebar_number_string);  
		// END OF Character Recognition Routine //}}}
	
	// Conversion of Units to Nanometres {{{
		if(scalebar_unit_string == "um"){
			scalebar_length_nm = scalebar_number_int * 1000;
		} else if(scalebar_unit_string == "nm"){
			scalebar_length_nm = scalebar_number_int;
		} // END OF Conversion of Units to Nanometres //}}}
	
	// Measure ScaleBar Length in Pixels //{{{
		w = getWidth();
		if(w == 1280){ M = 1; } else if(w == 1280*2){ M = 2; }
		msbl_y_line = 896*M+5;
		msbl_value_past = 0;
		msbl_tick_count = 0;
		for(x=0; x<w; x++){
			msbl_value_check = getPixel(x,msbl_y_line);
			if(msbl_value_check-msbl_value_past==255){
				msbl_last_switch_point = x;
				if(msbl_tick_count==0){
					msbl_first_switch_point = x;
				}
				msbl_tick_count++;
			}
			msbl_value_past = msbl_value_check;
		}
	
		scalebar_length_pix = msbl_last_switch_point -msbl_first_switch_point;
		pix_per_nm = scalebar_length_pix / scalebar_length_nm;
		nm_per_pixel = scalebar_length_nm / scalebar_length_pix;
		print(pix_per_nm + " px/nm \n" + nm_per_pixel + " nm/px"); //!@#$
		// END OF Measure ScaleBar Length in Pixels //}}}
		// END OF Measure Scale //}}}

		// Auto Crop //{{{
		run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000 frame=[0 sec] origin=0,0");
		w = getWidth(); if(w == 1280){ M = 1; } else if(w == 1280*2){ M = 2; } else { M = 0; }
		makeRectangle(0, 0, 1280*M, 896*M);
		run("Crop"); run("Select None"); //}}}
	}
	else exit("Something's wrong here");
	
	period_range_min_px = period_range_min_nm / nm_per_pixel; // For FFT period auto-detection
	period_range_max_px = period_range_max_nm / nm_per_pixel; // For FFT period auto-detection
	
	//OUTPUT DATA + TAGS
	w = getWidth(); h = getHeight();
	outputTD("nm_per_pixel",nm_per_pixel);
	outputTD("Width_initial",w);
	outputTD("Height_initial",h);

	
	// END OF DETERMINE RESOLUTION  //}}}
	
	// CROP STAGE 1 //{{{
	if((crop_choice == "Yes") & (crop_step == "Stage 1")){
		crop_choice_binary = 1;
		width_pc = getWidth();
		height_pc = getHeight();
		if(crop_uniform=="Uniform"){
			if(crop_unit == "Percent"){
				crop_x_left = crop_value / 100 * width_pc;
				crop_x_right = crop_value / 100 * width_pc;
				crop_y_top = crop_value / 100 * height_pc;
				crop_y_bottom = crop_value / 100 * height_pc;
			} else if(crop_unit == "Pixels") {
				crop_x_left = crop_value;
				crop_x_right = crop_value;
				crop_y_top = crop_value;
				crop_y_bottom = crop_value;
				
			}
		} else if(crop_uniform == "Non-Uniform"){
			if(crop_unit == "Percent"){
				crop_x_left = crop_x_left / 100 * width_pc;
				crop_x_right = crop_x_right / 100 * width_pc;
				crop_y_top = crop_y_top / 100 * height_pc;
				crop_y_bottom = crop_y_bottom / 100 * height_pc;
			}
		}
		// Ensure sufficient area: (Note, this could be modified to ensure area much larger than 1 pixel).
		if((crop_x_left+crop_x_right>=width_pc)||(crop_y_top+crop_y_bottom>=height_pc)){ exit("Your crop borders are too large.\n They don't leave anything!"); }
		// define crop box:
		xo_cbox = crop_x_left;                               
		yo_cbox = crop_y_top;
		w_cbox = width_pc - crop_x_left - crop_x_right;
		h_cbox = height_pc - crop_y_top - crop_y_bottom;
		makeRectangle(xo_cbox,yo_cbox,w_cbox,h_cbox);
		run("Duplicate...", "title="+"image001");
		image001 = getImageID();
		outputTD("Crop_1",crop_choice_binary);
		outputTD("Cxo",xo_cbox);
		outputTD("Cyo",yo_cbox);
		outputTD("Cw",w_cbox);
		outputTD("Ch",h_cbox);
	}
	else if(crop_choice == "No"){image001 = getImageID(); crop_choice_binary = 0; outputTD("Crop_1",crop_choice_binary);}
	selectImage(image001); saveAs("tiff",save_subfolder+"image001"+".tif");
	//}}}
	
	// PERIOD ESTIMATE BY WFFT  //{{{
	
	// A Macro to batch process images to measure spacings using FFT

	// Fast Fourier Transform
	run("FFT");
	img_fft = getImageID();
	
	// Crop Image to Reduce Time
	//    Time is proportionate to # of pixels!
	wfft = getWidth();
	downsizeto = round(0.25*wfft);
	if(w>=downsizeto){
	run("Canvas Size...", "width="+downsizeto+" height="+downsizeto+" position=Center zero");
	}
	
	
	// Do pixel-by-pixel intensity integration & counting
	
	Accuracy = 0.5;
	N = round(1/Accuracy);
	
	w = getWidth();
	h = getHeight();
	
	c = w/2;
	
	N_max = N*(sqrt(2)*c+1)+1;
	
	// Data Arrays
	ArrDist = newArray(N_max);
	ArrCnt = newArray(N_max);
	ArrInt = newArray(N_max);
	ArrAvgInt = newArray(N_max);
	
	for(n=0; n<N_max-1; n++){
	ArrDist[n] = n/N;
	ArrCnt[n] = 0;
	ArrInt[n] = 0;
	ArrAvgInt[n] = 0;
	}
	
	// Measures distance to centre for every pixel
	// Counts # of pixels at each distance
	// Integrates Intesnity as a Function of Distance
	for(y=0; y<h; y++){
	for(x=0; x<w; x++){
	d = sqrt(pow(x-c,2)+pow(y-c,2));
	val_a = getPixel(x,y);
	dd = round(d*N);
	val_o = ArrInt[dd];
	ArrInt[dd] = val_a+val_o;
	hits = ArrCnt[dd];
	hits++;
	ArrCnt[dd] = hits;
	}
	}
	
	// Calculates the Average Intensity
	// Also can fill-in the Results Table
	//run("Clear Results");
	for(n=0; n<N_max-1; n++){
	if(ArrCnt[n]>0){
	ArrAvgInt[n] = ArrInt[n] / ArrCnt[n];
	}
	//setResult("d",n,ArrDist[n]);
	//setResult("int",n,ArrInt[n]);
	//setResult("count",n,ArrCnt[n]);
	}
	
	
	//Piecewise smoothing Function
	for(n=0; n<N_max/2; n++){
	if(ArrAvgInt[n]==0){
	ArrAvgInt[n] = (ArrAvgInt[n+1]+ArrAvgInt[n-1])/2;
	}
	}
	for(n=N_max/2; n<N_max; n++){
	if(ArrAvgInt[n]==0){
	ArrAvgInt[n] = ArrAvgInt[n-1];
	}
	}

	// Smooth Values
	
	// Weighted Smooth Function
	function weightedsmooth(S,i,A,t) {
		if(i<A.length||i>S.length-A.length){return "ERROR";} else {
		initial = 0;
		for(n=0; n<A.length; n++){
			total = (S[i+n]+S[i-n])*A[n] + initial;
			initial = total;
		}
		return (total-A[0]*S[i])/t;
		}
	}
	
	// Weighted Smooth Normalizing Value
	function WSNV(A){
	initial = A[0];
	for(n=1; n<A.length; n++){
	total = initial + 2*A[n];
	initial = total;
	}
	return total;
	}
	
	weight_Array = newArray(8,7,6,5,4,3,2,1);
	wsnv = WSNV(weight_Array);
	
	ArrSmoothAvgInt = newArray(ArrAvgInt.length);
	for(n=0; n<ArrAvgInt.length; n++){
	ArrSmoothAvgInt[n] = ArrAvgInt[n];
	}
	
	for(n=weight_Array.length; n<ArrAvgInt.length-weight_Array.length; n++){
	ArrSmoothAvgInt[n] = weightedsmooth(ArrAvgInt,n,weight_Array,wsnv);
	}
  
	//Update Results Table
	wfft_real_distance = newArray(N_max);
	for(n=0; n<N_max-1; n++){
		wfft_real_distance[n] = wfft/ArrDist[n];
	}
	for(n=0; n<N_max-1; n++){
		setResult("real distance",n,wfft_real_distance[n]);
		setResult("radial distance",n,ArrDist[n]);
		setResult("int",n,ArrInt[n]);
		setResult("count",n,ArrCnt[n]);
		setResult("avg",n,ArrAvgInt[n]);
		setResult("smooth",n,ArrSmoothAvgInt[n]);
	}
	if(d_mode){updateResults(); selectWindow("Results"); saveAs("Measurements",save_subfolder+"d_mode_FFT.xls");}
  
	//   PEAK FINDING  
	n_peak = 0; intensity_peak = 0; d_peak_px = 0;
	for(n=0; n<N_max-1; n++){
		real_dist = wfft_real_distance[n];
		if((real_dist>period_range_min_px)&&(real_dist<period_range_max_px)){
			n_intensity = ArrSmoothAvgInt[n];
			if(n_intensity>intensity_peak){intensity_peak = n_intensity; n_peak = n; d_peak_px = real_dist;}
		}
	}
  
	// Calculate Final Values
	
	peak_position = ArrDist[n_peak];
	wfft_period_pixel = wfft / peak_position;
	wfft_period_nm = wfft_period_pixel * nm_per_pixel;
	print("WFFT Period: "+wfft_period_nm+" nm" + "\n" + "WFFT Period: "+wfft_period_pixel+" px");
	run("Clear Results");
	selectImage(img_fft); close();
	
	outputTD("wfft_Period_px",wfft_period_pixel);
	outputTD("wfft_Period_nm",wfft_period_nm);
	// END OF PERIOD ESTIMATE BY WFFT  //}}}
	
	
	
	// SMOOTHING //{{{
	//!@#$aul
	
	if(is("binary")){ outputTD("Smoothing",0);}
	else if(smoothing_method=="None"){ outputTD("Smoothing",0); }
	else{	
		if(smoothing_relative_option){smoothing_radius = smoothing_relative_percent/100*wfft_period_pixel;}
		run("Duplicate...", "title="+"image002");
		if(smoothing_method=="Automatic"){ outputTD("Smoothing",1);
			min_period_px = period_limits_nm[0]/nm_per_pixel; 
			max_period_px = period_limits_nm[1]/nm_per_pixel;
			auto_smoothing_radius = auto_smoothing_factor*wfft_period_pixel;
			//if(auto_smoothing_radius<=min_period_px || auto_smoothing_radius>=max_period_px){
			//	auto_smoothing_radius = auto_smoothing_factor*period_default_nm/nm_per_pixel;
			//}
			run("Median...", "radius="+auto_smoothing_radius); outputTD("Smoothing_radius",auto_smoothing_radius);
			if(d_mode){print("auto_smoothing_radius " + auto_smoothing_radius);}
			//!@#$ do I need to record this?
		}
		else if(smoothing_method=="Median"){ run("Median...", "radius="+smoothing_radius); outputTD("Smoothing",2); outputTD("Smoothing_radius",smoothing_radius);}
		else if(smoothing_method=="Gaussian"){ run("Gaussian Blur...", "sigma="+smoothing_radius); outputTD("Smoothing",3); outputTD("Smoothing_radius",smoothing_radius);}
	}
	image002 = getImageID(); saveAs("tiff",save_subfolder+"image002"+".tif"); 
	// END OF SMOOTHING //}}}

	// CROP STAGE 2 //{{{
	if((crop_choice == "Yes") & (crop_step == "Stage 2")){
		crop_choice_binary_two = 1;
		width_pc = getWidth();
		height_pc = getHeight();
		if(crop_uniform=="Uniform"){
			if(crop_unit == "Percent"){
				crop_x_left = crop_value / 100 * width_pc;
				crop_x_right = crop_value / 100 * width_pc;
				crop_y_top = crop_value / 100 * height_pc;
				crop_y_bottom = crop_value / 100 * height_pc;
			} else if(crop_unit == "Pixels") {
				crop_x_left = crop_value;
				crop_x_right = crop_value;
				crop_y_top = crop_value;
				crop_y_bottom = crop_value;
				
			}
		} else if(crop_uniform == "Non-Uniform"){
			if(crop_unit == "Percent"){
				crop_x_left = crop_x_left / 100 * width_pc;
				crop_x_right = crop_x_right / 100 * width_pc;
				crop_y_top = crop_y_top / 100 * height_pc;
				crop_y_bottom = crop_y_bottom / 100 * height_pc;
			}
		}
		// Ensure sufficient area: (Note, this could be modified to ensure area much larger than 1 pixel).
		if((crop_x_left+crop_x_right>=width_pc)||(crop_y_top+crop_y_bottom>=height_pc)){ exit("Your crop borders are too large.\n They don't leave anything!"); }
		// define crop box:
		xo_cbox = crop_x_left;                               
		yo_cbox = crop_y_top;
		w_cbox = width_pc - crop_x_left - crop_x_right;
		h_cbox = height_pc - crop_y_top - crop_y_bottom;
		makeRectangle(xo_cbox,yo_cbox,w_cbox,h_cbox);
		run("Duplicate...", "title="+"image003");
		image003 = getImageID();
		outputTD("Crop_2",crop_choice_binary_two);
		outputTD("Cxo",xo_cbox);
		outputTD("Cyo",yo_cbox);
		outputTD("Cw",w_cbox);
		outputTD("Ch",h_cbox);
	}
	else if(crop_choice == "No"){image003 = getImageID(); crop_choice_binary_two = 0; outputTD("Crop_2",crop_choice_binary_two);}
	selectImage(image003); saveAs("tiff",save_subfolder+"image003"+".tif");
	w = getWidth(); h = getHeight();
	//OUTPUT DATA + TAGS
	
	outputTD("Width_final",w);
	outputTD("Height_final",h);
	//}}}
	
	// THRESHOLDING //{{{
	run("Duplicate...", "title="+"image004");
	image004 = getImageID(); 
	if(is("binary")){ outputTD("Threshold",0); outputTD("Threshold.String","Binary");}
	else{
		outputTD("Threshold.String",thresholding_choice);
		//!@#$ Lots of stuff to implement here... could try my own bimodal histogram method
		if(thresholding_choice=="Custom"){ setAutoThreshold("Default dark"); output_threshold_type = 1; }
		else if(thresholding_choice=="Built-in"){
			setAutoThreshold("" + ij_thresholding_choice + " dark");
			outputTD("Threshold",2);
			outputTD("Threshold.Built-in.String",ij_thresholding_choice);
		}
		else if(thresholding_choice=="Auto-Local"){
		if(auto_local_relative){auto_local_radius = 1.5 * wfft_period_pixel;}
		run("Auto Local Threshold", "method="+auto_local_choice+" radius="+auto_local_radius+" parameter_1="+auto_local_p1+" parameter_2="+auto_local_p1+" white");
		outputTD("Threshold",3);
		outputTD("Threshold.Auto-Local.String",auto_local_choice);
		}
		
		
	}
	getThreshold(lower_threshold, upper_threshold); //!@#$
	outputTD("Up_Thresh",upper_threshold);
	outputTD("Low_Thresh",lower_threshold);
	print("Threshold\nLower: "+ lower_threshold + "\nUpper: " + upper_threshold);
	
	/* Following line, setThreshold(128,255), is added to change the behaviour 
	 * of the run("Make Binary") command. Without this line, the phases (bright and dark) 
	 * may be inverted, depending on the area fraction, such that the phase with the larger
	 * area fraction will always become the "positive" phase. 
	 * To swap the phases (i.e. to make the dark phase "positive"), the threshold can be
	 * changed using setThreshold(0,128) instead.
	 * cf: http://imagej.1557.x6.nabble.com/binary-image-inversion-irritation-td5004727.html
	 */
	setThreshold(128,255);
	run("Make Binary");
	saveAs("tiff",save_subfolder+"image004"+".tif"); 
	// END OF THRESHOLDING  //}}}

	// BINARY GROOMING //{{{
	// Functions to smooth the binary structures
	if(binary_grooming){
		
	edge_pts = w*2+h*2-4;
	xedge = newArray(edge_pts);
	yedge = newArray(edge_pts);
	i = 0;
	for(x=1; x<w; x++){ yedge[i] = 0; xedge[i] = x; i++; }
	for(y=1; y<h; y++){ yedge[i] = y; xedge[i] = w-1; i++; }
	for(x=w-2; x>=0; x--){yedge[i] = h-1; xedge[i] = x; i++; }
	for(y=h-2; y>=0; y--){yedge[i] = y; xedge[i] = 0; i++; }
	
	xi = newArray(round(w*h/4));
	yi = newArray(round(w*h/4));
	draw_erase = newArray(round(w*h/8));
	
	px_groom_condition = true;
	count = 0;
	ys = 1; ye = h-1;
	xs = 1; xe = w-1;
	
	
	
	next_time = false;
	now = false;
	
	while(px_groom_condition){
	//start = getTime();
		if(count>0){
			ymin = ye;
			ymax = ys;
			xmin = xe;
			xmax = ys;
		}
			
		for(i=0; i<count; i++){
			if(xi[i]>xmax){xmax = xi[i];}
			if(yi[i]>ymax){ymax = yi[i];}
			if(xi[i]<xmin){xmin = xi[i];}
			if(yi[i]<ymin){ymin = yi[i];}
			if(draw_erase[i]==1){setPixel(xi[i],yi[i],255);}
			if(draw_erase[i]==-1){setPixel(xi[i],yi[i],0);}
			draw_erase[i] = 0;
		}
		
		// Reframe:
		if(count>0){
			ys = maxOf(ymin-1,1);
			ye = minOf(ymax+1,h-1);
			xs = maxOf(xmin-1,1);
			xe = minOf(xmax+1,w-1);
		}
	
	i = 0;
	// Edges
	for(ei=0; ei<xedge.length; ei++){
		x = xedge[ei];
		y = yedge[ei];
		nv = neighbourValueExact(x,y,255);
		if(getPixel(x,y)==255){
			if(nv==0){
				draw_erase[i] = -1;
				xi[i] = x;
				yi[i] = y;
				i++;
			}
		} else if(nv>=3 && neighbourValueExactFour(x,y,255)>1){
			draw_erase[i] = 1;
			xi[i] = x;
			yi[i] = y;
			i++;
		}
	}
	
	for(y=ys; y<ye; y++){
		for(x=xs; x<xe; x++){
			if(getPixel(x,y)==0){ // DRAW
				nv = neighbourValueExact(x,y,255);
				//setPixel(x,y,nv);
				drawpx = false;
				if(nv>5){drawpx = true;}
				else if(nv>1) {
					cv = circuitValue(x,y,200);
					if(nv==2){ if(neighbourValueExactFour(x,y,255)>0&&cv==2){drawpx = true;} }
					else if(nv==5){
						if(neighbourValueExactFour(x,y,255)==3&&cv!=1){drawpx = true;}
						else if(cv==1&&neighbourValueExactFour(x,y,255)==3){ drawpx = true;}
					}
					else if(cv!=1){drawpx = true;}
				}
				if(drawpx){
					draw_erase[i] = 1;
					xi[i] = x;
					yi[i] = y;
					i++;
				}
			}
			else { // ERASE
				nv = neighbourValueExact(x,y,255);
				erasepx = false;
				if(nv<2){erasepx = true;}
				else if(nv<4){
					cv = circuitValue(x,y,200);
					if(nv==2){
						if(cv==1){erasepx = true;}
						else if(circuitValueX(x,y,200)==1 && circuitValueFour(x,y,200)==0){erasepx = true;}
					}
					else if(cv==1 && neighbourValueExactFour(x,y,255)==1){erasepx = true;}
				}
				if(erasepx){
					draw_erase[i] = -1;
					xi[i] = x;
					yi[i] = y;
					i++;
				}
			}
		}
		showProgress(y,h);
	}
	//time = getTime()-start; 
	//print(count + " " + time);
	if(i==count){px_groom_condition = false;}
	count = i;
	} // px_groom_condition LOOP
	
	// Edges AGAIN (Final) // NEEDS improvement and looping and orientation
	for(ei=0; ei<xedge.length; ei++){
		x = xedge[ei];
		y = yedge[ei];
		if(getPixel(x,y)==255){
			nv = neighbourValueExact(x,y,255);
			if(nv<=3 && neighbourValueExactFour(x,y,255)==1){setPixel(x,y,0);}
			else if(nv==0){setPixel(x,y,0);}
		}
	}
	
	}// End of "binary_grooming" //}}}

	// PARTICLE ANALYSIS //{{{
	// Cycle Loop //{{{
	cycles = 2;
	for(cyc=0; cyc<cycles; cyc++){
		cycle = cyc + 1;
		
		// Clean up the Binary Image //{{{
		/** IDEA: Should have different rules for thick vs thin images **/
		selectImage(image004);image_binary = getImageID(); //image_binary = image004;
		
		/** Noise elimination **/
		if(cyc==0){
			run("Set Measurements...", "area perimeter bounding shape redirect=None decimal=4");
			circ_limit = 0.35; //p.84 notebook
			small_area_limit = 0.75 * 3.14 * pow((wfft_period_pixel / 3),2);
			large_area_limit = 1.5 * 3.14 * pow((wfft_period_pixel / 2),2);
			for(i=0; i<2; i++){
				run("Invert");
				run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display clear record");
				setForegroundIndex(255);
				for(n=0;n<nResults;n++){
					if(getResult("Area",n)<small_area_limit){
						floodFill(getResult("XStart",n),getResult("YStart",n),"8-connected");
					}
					if(getResult("Area",n)<large_area_limit && getResult("Circ.",n)<circ_limit){
						floodFill(getResult("XStart",n),getResult("YStart",n),"8-connected");
					}
				}
			}
			
			/** Noise elimination -- Smoothing / mini pxMend **/ //!@#$ needs work
			for(i=0; i<2; i++){
				
				xmend = newArray(0); ymend = newArray(0);
				value = 255*(1-i);
				for(y=0; y<h; y++){
					for(x=0; x<w; x++){
						if(getPixel(x,y)==value){
							nv = neighbourValueExact(x,y,value);
							if(nv<3){ xmend = Array.concat(xmend,x); ymend = Array.concat(ymend,y);}
						}
					}
				}
				for(n=0; n<xmend.length; n++){ setPixel(xmend[n],ymend[n],255-value); }
			}
			imagesaver(savestages,image_binary,save_subfolder,"image004_pre","png");
		}
		// End of Clean-up//}}}
		
		// INITIAL Particle Analysis //{{{

		
		selectImage(image003);image_original = getImageID(); // image_original = image003;
		selectImage(image_original);
		run("Duplicate...", "title=positive");
		changeValues(0,0,1);
		image_positive = getImageID();
		selectImage(image_binary); run("Invert");
		imageCalculator("Subtract", image_positive,image_binary); 
		selectImage(image_binary); run("Invert");
		
		// Negative Image
		selectImage(image_original);
		run("Duplicate...", "title=negative");
		changeValues(255,255,254);
		image_negative = getImageID();
		imageCalculator("Add", image_negative,image_binary); 
		
							  
		//run("Set Measurements...", "area mean min perimeter bounding shape integrated redirect=None decimal=4");
		run("Set Measurements...", "area mean standard min perimeter bounding shape feret's integrated redirect=None decimal=4");
		
		
		// Analyze
		selectImage(image_positive);
		setThreshold(1,255);
		run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display clear record");
		nPositive = nResults();
		run("Create Selection");
		getRawStatistics(pos_nPixels, pos_mean, pos_min, pos_max, pos_std, pos_histogram);
		run("Select None");
			print("pos_nPixels " + pos_nPixels); 	outputTD("PA.pos_nPixels."+cycle,pos_nPixels);
			print("pos_mean " + pos_mean); 		outputTD("PA.pos_mean."+cycle,pos_mean);
			print("pos_min " + pos_min); 		outputTD("PA.pos_min."+cycle,pos_min);
			print("pos_max " + pos_max);		outputTD("PA.pos_max."+cycle,pos_max);
			print("pos_std " + pos_std);		outputTD("PA.pos_std."+cycle,pos_std);
		
		selectImage(image_negative);
		setThreshold(0,254);
		run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display record");
		nTotal = nResults();
		nNegative = nTotal - nPositive;
		run("Create Selection");
		getRawStatistics(neg_nPixels, neg_mean, neg_min, neg_max, neg_std, neg_histogram);
		run("Select None");
			print("neg_nPixels " + neg_nPixels); 	outputTD("PA.neg_nPixels."+cycle,neg_nPixels);
			print("neg_mean " + neg_mean); 		outputTD("PA.neg_mean."+cycle,neg_mean);
			print("neg_min " + neg_min); 		outputTD("PA.neg_min."+cycle,neg_min);
			print("neg_max " + neg_max);		outputTD("PA.neg_max."+cycle,neg_max);
			print("neg_std " + neg_std);		outputTD("PA.neg_std."+cycle,neg_std);
		
		selectImage(image_positive); imagesaver(savestages,image_positive,save_subfolder,"positive_"+cycle,"png"); close();
		selectImage(image_negative); imagesaver(savestages,image_positive,save_subfolder,"negative_"+cycle,"png"); close();
		
		print("nPositive = "+nPositive);	outputTD("PA.nPositive."+cycle,nPositive);
		print("nNegative = "+nNegative);	outputTD("PA.nNegative."+cycle,nNegative);
		print("nTotal = "+nTotal);		outputTD("PA.nTotal."+cycle,nTotal);
		
		
		// Mark phases
		for(n=0; n<nTotal; n++){
			if(n<nPositive){setResult("Phase",n,1);}
			else{setResult("Phase",n,0);}
		}
		
		// Check for Enclosure'
		selectImage(image_binary);
		checkInsideEDGE(nPositive,nTotal,100,w,h);
		
		// Edge Pixel Count
		selectImage(image_binary);
		edgePixelCount(100,0,0,w,h);
		
		//END of INITIAL //}}}
		
		// NOISE & Region-swapping
		/** IDEA: needs trimodal solution **/
		// Identifying regions with low/small differences between mean pixel intensities.
		if(cyc==0){
			selectImage(image_binary); //run("Duplicate...", "title=swapped"); image_swapped = getImageID();
			diff_threshold = minOf((pos_mean-neg_mean)/1.4,1.4*(pos_std+neg_std)); print("Diff_Threshold: " + diff_threshold); //!@#$
			setForegroundIndex(100);
			containing_table_entry = newArray(0);
			containing_table_entry_count = newArray(0);
			for(n=0; n<nResults; n++){
				if(getResult("Enclosed",n)==1){
					m = getResult("Enclosed.By",n);
					diff = abs(getResult("Mean",n)-getResult("Mean",m));
					setResult("Difference",n,diff);
					if(diff<diff_threshold){
						x = getResult("XStart",n); y = getResult("YStart",n);
						value = getPixel(x,y);
						setForegroundIndex(255-255*getResult("Phase",n)); floodFill(x,y,"8-connected");
						setResult("Swap",n,1); //phase = getResult("Phase",n); setResult("Phase",n,1-phase);
					}
		
				}
			}
		}
		updateResults();
		selectWindow("Results");
		saveAs("Measurements",save_subfolder+"feature_particle_analysis_"+cycle+".xls");
		imagesaver(savestages,image_binary,save_subfolder,"image004_cycle_"+cycle,"png");
	} // End of CYCLE Loop //}}}
		
		
		// WLSQ Perimeter & Area Fitting //{{{
		//#!# needs to be conditional on the number of results
		large_drops = 2; //#!# user setting???
		min_area = wfft_period_pixel*wfft_period_pixel*0.125; //#!# user setting??? assumes small = 1/4 of period
		wlsq_iterations = 5; //#!# user setting???
			outputTD("PA.SET.large_drops",large_drops);
			outputTD("PA.SET.min_area",min_area);
			outputTD("PA.SET.wlsq_iterations",wlsq_iterations);
		selectImage(image_binary);
		positive_wlsq = particleWLSQfull(0,nPositive,"median_distance",wlsq_iterations,min_area,large_drops);
		negative_wlsq = particleWLSQfull(nPositive,nTotal,"median_distance",wlsq_iterations,min_area,large_drops);
		Array.print(positive_wlsq);
		Array.print(negative_wlsq);
		outputTD("PA.WLSQ.positive.0",positive_wlsq[0]);
		outputTD("PA.WLSQ.positive.1",positive_wlsq[1]);
		outputTD("PA.WLSQ.positive.2",positive_wlsq[2]);
		outputTD("PA.WLSQ.positive.3",positive_wlsq[3]);
		outputTD("PA.WLSQ.positive.4",positive_wlsq[4]);
		outputTD("PA.WLSQ.positive.5",positive_wlsq[5]);
		outputTD("PA.WLSQ.positive.6",positive_wlsq[6]);
		outputTD("PA.WLSQ.negative.0",negative_wlsq[0]);
		outputTD("PA.WLSQ.negative.1",negative_wlsq[1]);
		outputTD("PA.WLSQ.negative.2",negative_wlsq[2]);
		outputTD("PA.WLSQ.negative.3",negative_wlsq[3]);
		outputTD("PA.WLSQ.negative.4",negative_wlsq[4]);
		outputTD("PA.WLSQ.negative.5",negative_wlsq[5]);
		outputTD("PA.WLSQ.negative.6",negative_wlsq[6]);
		
		// Generating usable, failsafe values for widths
		width_accuracy_est = 0.985;
		if((positive_wlsq[0] == -1) & (negative_wlsq[0] == -1)){
			// If both fail (something's horribly wrong, but what the heck... assume 50-50)
			negative_width = wfft_period_pixel/2;
			positive_width = wfft_period_pixel/2;
			proportion_width = 0.5;
		}
		else if((positive_wlsq[0] != -1) & (negative_wlsq[0] == -1)){ 
			// If positive fails (i.e. only 1 piece / black dots)
			positive_width = positive_wlsq[0]*width_accuracy_est;
			negative_width = wfft_period_pixel-positive_width;
			proportion_width = positive_width/(negative_width+positive_width);
		}
		else if((positive_wlsq[0] == -1) & (negative_wlsq[0] != -1)){ 
			// If negative fails (i.e. only 1 piece / white dots)  
			negative_width = negative_wlsq[0]*width_accuracy_est;
			positive_width = wfft_period_pixel-negative_width;
			proportion_width = positive_width/(negative_width+positive_width);
		}
		else {	// if everything works:
			proportion_width = positive_wlsq[0]/(negative_wlsq[0]+positive_wlsq[0]);
			positive_width = wfft_period_pixel*proportion_width;
			negative_width = wfft_period_pixel*(1-proportion_width);
		}
		outputTD("PA.Width.positive",positive_width); print("PA.Width.positive: "+positive_width); 
		outputTD("PA.Width.negative",negative_width); print("PA.Width.negative: "+negative_width);
		outputTD("PA.Width.proportion",proportion_width); print("PA.Width.proportion: "+proportion_width);
		// End of WLSQ //}}}
		
		
		
		// FEATURE Analysis //{{{
		
		// BLOBS//{{{
		run_blobs = false;
		if(run_blobs){
		// Positive Blobs
		B_p_thresh = 1.5 * positive_width; pB_area_min = 4*positive_width*positive_width;
		selectImage(image_binary);
		run("Duplicate...", "title=blobs_pos");
		run("Distance Map");
			changeValues(B_p_thresh,255,255);changeValues(0,B_p_thresh,0);
		image_blobs_pos = getImageID(); imagesaver(savestages,image_blobs_pos,save_subfolder,"image_blobs_pos_i","png");
		
		pB_area = newArray(0); pB_x = newArray(0); pB_y = newArray(0); area_sum = 0;
		for(y=0; y<h; y++){ for(x=0; x<w; x++){
			if(getPixel(x,y)==255){
				doWand(x,y); getStatistics(area); changeValues(255,255,100);
				run("Select None"); area_sum += area;
				pB_x = Array.concat(pB_x,x); pB_y = Array.concat(pB_y,y); pB_area = Array.concat(pB_area,area);
			}
		}}
		changeValues(100,100,255);
		outputTD("PA.Blobs.p.i.count",pB_area.length); outputTD("PA.Blobs.p.i.area",area_sum);
		setForegroundIndex(0); Blobs_erased = 0; 
		for(i=0; i<pB_area.length; i++){ if(pB_area[i]< pB_area_min){floodFill(pB_x[i],pB_y[i],"8-connected"); Blobs_erased++; area_sum -= pB_area[i];}}
		pos_blob_count = pB_area.length-Blobs_erased;  outputTD("PA.Blobs.p.f.area",area_sum);
		outputTD("PA.Blobs.p.f.count",pos_blob_count); imagesaver(savestages,image_blobs_pos,save_subfolder,"image_blobs_pos_f","png");
		
		// Negative Blobs
		B_n_thresh = 1.5 * negative_width; nB_area_min = 4*negative_width*negative_width;
		selectImage(image_binary);
		run("Duplicate...", "title=blobs_neg");
		run("Invert");
		run("Distance Map");
			changeValues(B_n_thresh,255,255);changeValues(0,B_n_thresh,0);
		image_blobs_neg = getImageID(); imagesaver(savestages,image_blobs_neg,save_subfolder,"image_blobs_pos_i","png");
		
		nB_area = newArray(0); nB_x = newArray(0); nB_y = newArray(0); area_sum = 0;
		for(y=0; y<h; y++){ for(x=0; x<w; x++){
			if(getPixel(x,y)==255){
				doWand(x,y); getStatistics(area); changeValues(255,255,100);
				run("Select None"); area_sum += area;
				nB_x = Array.concat(nB_x,x); nB_y = Array.concat(nB_y,y); nB_area = Array.concat(nB_area,area);
			}
		}}
		changeValues(100,100,255);
		outputTD("PA.Blobs.n.i.count",nB_area.length); outputTD("PA.Blobs.n.i.area",area_sum);
		setForegroundIndex(0); Blobs_erased = 0; 
		for(i=0; i<nB_area.length; i++){ if(nB_area[i]< nB_area_min){floodFill(nB_x[i],nB_y[i],"8-connected"); Blobs_erased++; area_sum -= nB_area[i];}}
		neg_blob_count = nB_area.length-Blobs_erased;  outputTD("PA.Blobs.n.f.area",area_sum);
		outputTD("PA.Blobs.n.f.count",neg_blob_count); imagesaver(savestages,image_blobs_neg,save_subfolder,"image_blobs_pos_f","png");
		}//}}} BLOBS
		
		// DOTS
		/** IDEA: Needs a more nuanced approach for the dots **/
		dot_max_area_pos = 6*PI*pow(positive_width/2,2);
		pos_edge_dot_count = 0; pos_dot_count = 0;
		dot_max_area_neg = 6*PI*pow(negative_width/2,2);
		neg_edge_dot_count = 0; neg_dot_count = 0;
		selectImage(image_binary); run("Duplicate...", "title=dots_pos"); image_dots_pos = getImageID();
		selectImage(image_binary); run("Duplicate...", "title=dots_neg"); image_dots_neg = getImageID(); run("Invert");
		setForegroundIndex(0);
		for(n=0; n<nResults; n++){
			if(getResult("Phase",n)==1){
				if(getResult("Area",n) <= dot_max_area_pos){
					setResult("Dot",n,1);
					if(getResult("OnEdge",n)==1){pos_edge_dot_count += 1;} else {pos_dot_count += 1;}
				} else {
					setResult("Dot",n,0);
					selectImage(image_dots_pos);floodFill(getResult("XStart",n),getResult("YStart",n),"8-connected");
				}
			}
				if(getResult("Phase",n)==0){	
					if(getResult("Area",n) <= dot_max_area_neg){
						setResult("Dot",n,1);
						if(getResult("OnEdge",n)==1){neg_edge_dot_count += 1;} else {neg_dot_count += 1;}
				} else {
					setResult("Dot",n,0);
					selectImage(image_dots_neg);floodFill(getResult("XStart",n),getResult("YStart",n),"8-connected");
				}
			}
		}
		outputTD("dot_max_area_pos",dot_max_area_pos);
		outputTD("pos_edge_dot_count",pos_edge_dot_count);
		outputTD("pos_dot_count",pos_dot_count);
		outputTD("dot_max_area_neg",dot_max_area_neg);
		outputTD("neg_edge_dot_count",neg_edge_dot_count);
		outputTD("neg_dot_count",neg_dot_count);
		
		// DOT COORDINATES
		raw_dot_count = pos_edge_dot_count + pos_dot_count + neg_edge_dot_count + neg_dot_count;
		//dots_xy = newArray(raw_dot_count); 
		dot_defects = newArray(raw_dot_count);
		//dot_orientations = newArray(raw_dot_count); //Not implemented yet; placeholder
		di = 0;
		for(n=0; n<nResults; n++){
			if( (getResult("Dot",n)==1) && (getResult("OnEdge",n)==0) ){
				phase = getResult("Phase");
				x = round(getResult("BX",n) + getResult("Width",n)/2); //change particle analysis to make it the centre of mass?
				y = round(getResult("BY",n) + getResult("Height",n)/2);
				dot_defects[di] = defect_coder(0,phase,0,x,y,"enc");
				di++;
			}
		}
		for(n=0; n<nResults; n++){
			if( (getResult("Dot",n)==1) && (getResult("OnEdge",n)==1) ){
				phase = getResult("Phase");
				x = round(getResult("BX",n) + getResult("Width",n)/2); //change particle analysis to make it the centre of mass?
				y = round(getResult("BY",n) + getResult("Height",n)/2);
				dot_defects[di] = defect_coder(0,phase,1,x,y,"enc");
				di++;
			}
		}
		
		// LINES
		/** IDEA: Needs a more nuanced approach for the lines **/
		line_min_area_pos = 2.5*PI*pow(positive_width/2,2);
		pos_edge_line_count = 0; pos_line_count = 0;
		line_min_area_neg = 2.5*PI*pow(negative_width/2,2);
		neg_edge_line_count = 0; neg_line_count = 0;
		selectImage(image_binary); run("Duplicate...", "title=lines_pos"); image_lines_pos = getImageID();
		selectImage(image_binary); run("Duplicate...", "title=lines_neg"); image_lines_neg = getImageID(); run("Invert");
		setForegroundIndex(0);
		for(n=0; n<nResults; n++){
			if(getResult("Phase",n)==1){
				if(getResult("Area",n) > line_min_area_pos){
					setResult("Line",n,1);
					if(getResult("OnEdge",n)==1){pos_edge_line_count += 1;} else {pos_line_count += 1;}
				} else {
					setResult("Line",n,0);
					selectImage(image_lines_pos);floodFill(getResult("XStart",n),getResult("YStart",n),"8-connected");
				}
			}
				if(getResult("Phase",n)==0){	
					if(getResult("Area",n) > line_min_area_neg){
						setResult("Line",n,1);
						if(getResult("OnEdge",n)==1){neg_edge_line_count += 1;} else {neg_line_count += 1;}
				} else {
					setResult("Line",n,0);
					selectImage(image_lines_neg);floodFill(getResult("XStart",n),getResult("YStart",n),"8-connected");
				}
			}
		}
		outputTD("line_min_area_pos",line_min_area_pos);
		outputTD("pos_edge_line_count",pos_edge_line_count);
		outputTD("pos_line_count",pos_line_count);
		outputTD("line_min_area_neg",line_min_area_neg);
		outputTD("neg_edge_line_count",neg_edge_line_count);
		outputTD("neg_line_count",neg_line_count);
		

		// Identify dots patches

		
		selectWindow("Results");
		saveAs("Measurements",save_subfolder+"feature_particle_analysis.xls");
		//run("Clear Results");
		//}}}

		// Mark Phases  //!@#$

		// Repeat Period Calculation  //!@#$
	
	
	// END OF PARTICLE ANALYSIS //}}}
	
	
	// DIVIDE IMAGES //{{{
	/* Need to consider what to do if no particles for one or other */
	pos_min_area = line_min_area_pos;
	neg_min_area = line_min_area_neg;
	image004 = image_binary;
	
		// POSITIVE LINES
		selectImage(image004); //binary //#!# Temporary
		run("Duplicate...", "title=image_positive_lines");
		image_positive_lines = getImageID();
		colourParticlesPhase(1,"Area","<",pos_min_area,0,0,nResults);
		if(run_blobs){imageCalculator("Subtract", image_positive_lines,image_blobs_pos);}
		
		// POSITIVE DOTS
		selectImage(image004); //binary //#!# Temporary
		run("Duplicate...", "title=image_positive_dots");
		image_positive_dots = getImageID();
		colourParticlesPhase(1,"Area",">=",pos_min_area,0,0,nResults);
		
		// NEGATIVE LINES
		selectImage(image004); //binary //#!# Temporary
		run("Duplicate...", "title=image_negative_lines");
		image_negative_lines = getImageID();
		run("Invert");
		colourParticlesPhase(0,"Area","<",neg_min_area,0,0,nResults);
		if(run_blobs){imageCalculator("Subtract", image_negative_lines,image_blobs_neg);}
		
		// NEGATIVE DOTS
		selectImage(image004); //binary //#!# Temporary
		run("Duplicate...", "title=image_negative_dots");
		image_negative_dots = getImageID();
		run("Invert");
		colourParticlesPhase(0,"Area",">=",neg_min_area,0,0,nResults);
		outputTD("pos_min_area",pos_min_area); outputTD("neg_min_area",neg_min_area);
		
	
	// END of DIVIDE IMAGES //}}}
	
	// SKELETON ANALYSES // {{{
		
		// CREATE SKELETONS //{{{
		//Positive Skeleton
		/**run("Make Binary"); // for good measure?**/
		selectImage(image_positive_lines);
		run("Duplicate...", "title=image_positive_skeleton");
		run("Skeletonize");
		image_positive_skeleton = getImageID();
		
		// NEGATIVE SKELETON
		/**run("Make Binary"); // for good measure?**/
		selectImage(image_negative_lines);
		run("Duplicate...", "title=image_negative_skeleton");
		run("Skeletonize");
		image_negative_skeleton = getImageID();
		
		// END of CREATE SKELETONS //}}}
		
		// PHASE ONE
		
		// SKELETON GROOMING //{{{
		
		for(i_phase = 0; i_phase < 2; i_phase++){
			if(i_phase==0){
				selectImage(image_positive_skeleton);
				mDist = positive_width; outputTD("pos_mDist",mDist);
			}
			else{
				selectImage(image_negative_skeleton);
				mDist = negative_width; outputTD("neg_mDist",mDist);
			}
		
// Grooming_New_006.txt
  //itime = getTime();
  //mDist = Width_Est_r;
  
  Li_array = newArray(-1,-1,-1,0,1,1,1,0);
  Lj_array = newArray(-1,0,1,1,1,0,-1,-1);
  
  w = getWidth();
  h = getHeight();

  // FIND TERMINAL POINTS & JUNCTIONS [combined search]
  showStatus("Groom: Finding Junctions & Terminal Points");

  n = 0;
  m = 0;
  for (y=0; y<h; y++) {
    for (x=0; x<w; x++) {
      if (getPixel(x,y)==255){
        ptvalue = getPixel(x-1,y-1)+getPixel(x-1,y)+getPixel(x-1,y+1)+getPixel(x,y+1)+getPixel(x+1,y+1)+getPixel(x+1,y)+getPixel(x+1,y-1)+getPixel(x,y-1);
        if(ptvalue==255) {
          setResult("X_TP",n,x);
          setResult("Y_TP",n,y);
          n++;
        }
        if(ptvalue>2*255) {
          setResult("X_JP",m,x);
          setResult("Y_JP",m,y);
          m++;
        }
      }
    }
    if (y%10==0) showProgress(y,h);
  }
  
  nTP = n;
  TPval_X = newArray(nTP);
  TPval_Y = newArray(nTP);
  for(n=0; n<nTP; n++){
    TPval_X[n] = getResult("X_TP",n);
    TPval_Y[n] = getResult("Y_TP",n);
  }
  TP_marked = newArray(nTP);
  
  nJP_pre = m;
  JPval_X_pre = newArray(nJP_pre);
  JPval_Y_pre = newArray(nJP_pre);
  for(n=0; n<nJP_pre; n++){
    JPval_X_pre[n] = getResult("X_JP",n);
    JPval_Y_pre[n] = getResult("Y_JP",n);
    setPixel(JPval_X_pre[n],JPval_Y_pre[n],10);   // change from 10 to 0
  }


  run("Clear Results");
    
  // START: Remove single-pixel tags.//{{{
  
  // This section checks for 2-connected pixels which are attached to Junction points, in isolation
  // but which are not themselves part of the junction
  showStatus("Groom: Trim Single-Pixel Tags");
  tagMend = true;
  if(tagMend){
    Jval = 10;
    wr = getWidth()-1;
    hr = getHeight()-1;
  
    for (y=1; y<hr; y++) {
        for (x=1; x<wr; x++) {
        if (getPixel(x,y)==255){
            ptvalue = floor((getPixel(x-1,y-1)+getPixel(x-1,y)+getPixel(x-1,y+1)+getPixel(x,y+1)+getPixel(x+1,y+1)+getPixel(x+1,y)+getPixel(x+1,y-1)+getPixel(x,y-1))/Jval);
                if(ptvalue==2) {
                    ptvalue_x = floor((getPixel(x-1,y)+getPixel(x+1,y))/Jval);
                    ptvalue_y = floor((getPixel(x,y-1)+getPixel(x,y+1))/Jval);
                    ptvalue_top = floor((getPixel(x-1,y+1)+getPixel(x,y+1)+getPixel(x+1,y+1))/Jval);
                    ptvalue_left = floor((getPixel(x-1,y+1)+getPixel(x-1,y)+getPixel(x-1,y-1))/Jval);
                    switch1 = 0;
                    switch2 = 0;
                    
                    if(ptvalue_y==0||ptvalue_x==0){
  
                        // X-axis check
  
                        if(ptvalue_x==0){
                            if(ptvalue-ptvalue_top>=1&&ptvalue_top!=0){
                            switch1++;
                            }
                            else{
                            switch2++;
                            }
  
                        }
  
                        // Y-axis check
  
                        if(ptvalue_y==0){
                            if(ptvalue-ptvalue_left>=1&&ptvalue_left!=0){
                            switch1++;
                            }
                            else{
                            switch2++;
                            }
                        }
  
                    }
  
                    else {
                    ptvalue_L = getPixel(x+1,y)+getPixel(x,y+1);
                    ptvalue_F = getPixel(x+1,y)+getPixel(x,y-1);
                    ptvalue_Z = getPixel(x-1,y)+getPixel(x,y-1);
                    ptvalue_J = getPixel(x-1,y)+getPixel(x,y+1);
  
                        if(ptvalue_L<Jval){
                            if(getPixel(x+1,y+1)==Jval){
                            switch1++;
                            }
                        }
                        if(ptvalue_F<Jval){
                            if(getPixel(x+1,y-1)==Jval){
                            switch1++;
                            }
                        }
                        if(ptvalue_Z<Jval){
                            if(getPixel(x-1,y-1)==Jval){
                            switch1++;
                            }
                        }
                        if(ptvalue_J<Jval){
                            if(getPixel(x-1,y+1)==Jval){
                            switch1++;
                            }
                        }
  
                    }
  
  
                        if(switch1-switch2>=1){
                        setPixel(x,y,25);
                        }
                        else { setPixel(x,y,5);}
  
                    
                }
            }
        }
        if (y%10==0) showProgress(y,hr);
    }
  }
  
  changeValues(5,5,0);
  changeValues(10,25,255);
  run("Skeletonize");
  
  // END: Remove single-pixel tags. //}}}
  
  // Re-define junctions:
  m = 0;
  for(n=0; n<nJP_pre; n++){
    x = JPval_X_pre[n];
    y = JPval_Y_pre[n];
    if(getPixel(x,y)==255){
      m++;
    }
  }
  nJP = m;
  
  JPval_X = newArray(nJP);
  JPval_Y = newArray(nJP);
  
  m = 0;
  for(n=0; n<nJP_pre; n++){
    x = JPval_X_pre[n];
    y = JPval_Y_pre[n];
    if(getPixel(x,y)==255){
      JPval_X[m] = x;
      JPval_Y[m] = y;
      m++;
    }
  }

  JP_TPcount = newArray(nJP);
  JP_nearestTP = newArray(nJP);
  JP_Deg_G = newArray(nJP);
  
  
  // ####
  
  // JUNCTION DEGREE
  // Determine Junction Connectivity
  showStatus("Groom: Junction Types");
  MaxPx = 255;
  for(m=0; m<nJP; m++){
    x = JPval_X[m];
    y = JPval_Y[m];
    ptvalue = getPixel(x-1,y-1)+getPixel(x-1,y)+getPixel(x-1,y+1)+getPixel(x,y+1)+getPixel(x+1,y+1)+getPixel(x+1,y)+getPixel(x+1,y-1)+getPixel(x,y-1);
    JP_Deg_G[m] = ptvalue/MaxPx;
  }
  // Colouring the Junctions
  for(m=0; m<nJP; m++){
    setPixel(JPval_X[m],JPval_Y[m],JP_Deg_G[m]+1);
  }
  
  // Calculating Junction Degree
  run("Set Measurements...", "area mean center bounding display redirect=None decimal=4");
  setThreshold(1, 5);
  run("Analyze Particles...", "size=1-Infinity circularity=0.00-1.00 show=Nothing display");
  resetThreshold();
  nJ_pts = nResults();
  changeValues(1,10,10);  // changed back to 10
  Jn_X = newArray(nJ_pts);
  Jn_Y = newArray(nJ_pts);
  Jn_D_extra = newArray(nJ_pts);
  for(n=0; n<nJ_pts; n++){
    Jn_X[n] = getResult("XM",n);
    Jn_Y[n] = getResult("YM",n);
    jArea = getResult("Area",n);
    jMean = getResult("Mean",n);
    Jn_D_extra[n] =round(jArea*(jMean-1))-2;
  }
  run("Clear Results");
  
  
  // Find the nearest TP & Mark it.
  showStatus("Groom: Nearest Terminal Points");
  for (k=0; k<nJ_pts; k++) {
    i = Jn_X[k];
    j = Jn_Y[k]; 
    minDist = mDist+2;
    n_coterminals = 0;
    KT = 0;
    for (t=0; t<nTP; t++) {
      x = TPval_X[t];
      y = TPval_Y[t];
      Dist = sqrt(pow(i-x,2)+pow(j-y,2));
      if(Dist<=mDist+1){
          n_coterminals++;
          if(Dist<minDist){
            minDist = Dist;
            KT = t;
          }
        }
      }
       if(n_coterminals>0){
      TP_marked[KT] = 1;
       }
    }
 
  // ####
  
  // MEASURE JUNCTION-TERMINAL DISTANCES
  // Formerly finding "coterminal points"
  //... remove short branches
  showStatus("Groom: Measuring JP-TP Distances");
  
  coterminal_Js = 0;
  for (k=0; k<nJP; k++) {
    i = JPval_X[k];
    j = JPval_Y[k]; 
    minDist = mDist+1;
    n_coterminals = 0;
    KT = 0;
    for (t=0; t<nTP; t++) {
      mark = TP_marked[t];
      if(mark==1){
        x = TPval_X[t];
        y = TPval_Y[t];
        Dist = sqrt(pow(i-x,2)+pow(j-y,2));
        if(Dist<=mDist){
          pijvalue1 = getPixel(i-1,j-1)+getPixel(i-1,j)+getPixel(i-1,j+1)+getPixel(i,j+1)+getPixel(i+1,j+1)+getPixel(i+1,j)+getPixel(i+1,j-1)+getPixel(i,j-1);
          setColor(55);
          floodFill(x,y,"8-connected");
          pijvalue2 = getPixel(i-1,j-1)+getPixel(i-1,j)+getPixel(i-1,j+1)+getPixel(i,j+1)+getPixel(i+1,j+1)+getPixel(i+1,j)+getPixel(i+1,j-1)+getPixel(i,j-1);
          if(pijvalue1-pijvalue2==200){
            n_coterminals++;
            if(Dist<minDist){
              minDist = Dist;
              KT = t;
            }
          }
          setColor(255);
          floodFill(x,y,"8-connected");  
        }
      }
    }
    if(n_coterminals>0){
      coterminal_Js++;
    }
    JP_nearestTP[k] = KT;
    JP_TPcount[k] = n_coterminals;
  }
  
  for (k=0; k<nJP; k++) {
    i = JPval_X[k];
    j = JPval_Y[k];
    if(getPixel(i,j)==10){
      if(getPixel(i-1,j-1)>=10){ ijA=1;} else {ijA =0;}
      if(getPixel(i,j-1)>=10){ijB=1;} else {ijB =0;}
      if(getPixel(i+1,j-1)>=10){ijC=1;} else {ijC =0;}
      if(getPixel(i+1,j)>=10){ijD=1;} else {ijD =0;}
      if(getPixel(i+1,j+1)>=10){	ijE=1;} else {ijE =0;}
      if(getPixel(i,j+1)>=10){	ijF=1;} else {ijF =0;}
     if(getPixel(i-1,j+1)>=10){ijG=1;} else {ijG =0;}
      if(getPixel(i-1,j)>=10){ijH=1;} else {ijH =0;}
      pijNvalue = ijA+ijB+ijC+ijD+ijE+ijF+ijG+ijH;
      if(pijNvalue>=3){
        if(JP_TPcount[k]>0){
          a = JP_nearestTP[k];
          x = TPval_X[a];
          y = TPval_Y[a];
          if(getPixel(x,y)==255){
            
            pijvalue1 = getPixel(i-1,j-1)+getPixel(i-1,j)+getPixel(i-1,j+1)+getPixel(i,j+1)+getPixel(i+1,j+1)+getPixel(i+1,j)+getPixel(i+1,j-1)+getPixel(i,j-1);
            setColor(55);
            floodFill(x,y,"8-connected");
            pijvalue2 = getPixel(i-1,j-1)+getPixel(i-1,j)+getPixel(i-1,j+1)+getPixel(i,j+1)+getPixel(i+1,j+1)+getPixel(i+1,j)+getPixel(i+1,j-1)+getPixel(i,j-1);
            if(pijvalue1-pijvalue2==200){
              setColor(0);
              floodFill(x,y,"8-connected");
              setColor(20);
              floodFill(i,j,"8-connected");
              // Now to determine whether to delete ij
              pijvalue3 = getPixel(i-1,j-1)+getPixel(i-1,j)+getPixel(i-1,j+1)+getPixel(i,j+1)+getPixel(i+1,j+1)+getPixel(i+1,j)+getPixel(i+1,j-1)+getPixel(i,j-1);
              if(pijvalue3<255){
                setPixel(i,j,0);
              }
            }
            else{
              setColor(255);
              floodFill(x,y,"8-connected");
            }
            A = JP_TPcount[k];
            A--;
            JP_TPcount[k] = A;
          }
        }
      }
    }
  }
  
  changeValues(10,20,255);  
  run("Skeletonize");
  //ftime = getTime();
  //print((ftime-itime)/1000);

    showStatus("Groom: Complete!");		
		
		} // End of i_phase loop
		
/** Previous versions had the following sections:

Before skeletonization:
// Border Pixel Reconnect
// Bringing Terminal Points near the Edge to the Edge
// Exclusive @BCPMM

// ##### SG 01 #####
// ELIMINATE HOLES
// Cleaning up small holes in the skeleton

// ##### SG 02 #####
// GROOM SKELETON  
// Grooming_New_006.txt

// ##### SG 03 ##### 
// PERIOD MEASUREMENT
	// REMOVE JUNCTIOS & LOOPS
	// ## START: _Junction_and_Loop_Remover_001.txt
	// PERIOD MEASUREMENT
	// ## START: Wire_Periodicity_007.txt
	
	Extracted this as: SG_03_Wire_Period_Calc_extraced_20140115.txt
	Works fast for small areas; scales exponentially for full image. 
	This could be a quick back-up method for confirming period 
	
	Would need to limit area to parts well covered by lines.
	Centre... if coverage > 90%, anywhere works.
	Could also use in small patches over the entire image. 8*12 sec <<< 70 min
	
	Contains a useful element for finding loops w/o junctions.
	> could extend to finding nested loops

// ##### SG 04 #####  (erased this step in 100e_38_express)
// EXTEND TERMINAL POINTS

// ##### SG 05 ##### (Optional in 100e_38_express)
// RECONNECTION
**/
		
		// END OF SKELETON GROOMING //}}}
		
		// CHARACTERIZE JUNCTIONS & TERMINAL POINTS //{{{
	
	defects_array = newArray(0); // This will contain all defects, valid and non-valid.
	defects_orientations = newArray(0); // This will contain all defect orientations
	// First Positive, Then Negative
	
	for(i_phase = 0; i_phase < 2; i_phase++){
		if(i_phase==0){
			selectImage(image_positive_skeleton);
			run("Duplicate...", "title=positive_j_skeleton");
			positive_j_skeleton = getImageID();
			phase = 0;
		}
		else{
			selectImage(image_negative_skeleton);
			run("Duplicate...", "title=negative_j_skeleton");
			negative_j_skeleton = getImageID();
			phase = 1;
		}
	
	// Junction Function (What's Your Conjunction?)
	// Using doWand, getStatistics, & getSelectionCoordinates

	// Step 1: Mark the skeleton
	// Marking Guide: (2)-connected -> 2; (3...8)-conected -> 10

	// First sweep: Initial Marks
	/** 	Marking Guide: 
		(2)-connected -> 2;
		These will need to be subsequently marked. (Second sweep)
		(3...8)-conected -> 10
		**/
	
	w = getWidth(); h = getHeight();
	
	t_count = 0;
	for(y=0; y<h; y++){
		for(x=0; x<w; x++){
			if(getPixel(x,y)==255){
				n_value = neighbourValue(x,y,1);
				if(n_value >= 3){ n_value = 10; } // ensures junction pts all 3
				setPixel(x,y,n_value);
				if(n_value == 1){ t_count++; }
			}
		}
		showProgress(y,h);
	}
	
	// Second sweep: Specific Marks
	/** 	2-connected adjacent to 1 junction piyel -> 9
		2-connected adjacent to 2 junction piyels -> 11
		Note that a 2-connected piece cannot be
		adjacent to more than 2 junction piyels by definition.
		The second case is rare but possible.
		Combined this enables analysis of the junction **/
	
	terminal_xy = newArray(t_count); ti = 0;
	terminal_defects = newArray(t_count);
	terminal_orientations = newArray(t_count); //Not implemented yet; placeholder
	for(y=0; y<h; y++){
		for(x=0; x<w; x++){
			if(getPixel(x,y)==2){
				n_value = neighbourValueExact(x,y,10);
				if(n_value == 1){ 
					px_value = 9; 
					setPixel(x,y,px_value);
				}
				if(n_value > 1){ 
					px_value = 11; 
					setPixel(x,y,px_value);
				}
				
			}
			if(getPixel(x,y)==1){
				terminal_xy[ti] = xy_coder(x,y,0,"enc");
				terminal_defects[ti] = defect_coder(0,phase,1,x,y,"enc");
				terminal_orientations[ti] = 0; //TEMPORARY
				ti++; 
			}
		}
		showProgress(y,h);
	}
	
	// Step 2: Measure the Junctions
	
	j_count = 0;
	junctions = newArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	start = getTime();
	junction_xy = newArray(0);
	junction_type = newArray(0);
	junction_defects = newArray(0);
	junction_orientations = newArray(0); //v0.50d new addition
	changeValues(9,9,12); // v029... prevents bridging
	//w = getWidth(); h = getHeight(); 
	for(y=0; y<h; y++){
		for(x=0; x<w; x++){
			if(getPixel(x,y)==10){
				doWand(x,y,1,"8-connected"); // tolerance set to 1: 3
				//getStatistics(area, mean, min, max, std, histogram);
				j_count++;
				//jtype = histogram[9];
				jtype = junctionDegree(12,13);
				junctions[jtype] += 1;
				//getSelectionCoordinates(xpoints,ypoints);
				//xy = centrePixel(xpoints,ypoints,1);
				// Better-centred:
				//changeValues(9,9,12);
				//doWand(x,y,1,"8-connected");
				xy = getSelectionPixels(3);
				junction_xy = Array.concat(junction_xy,xy);
				junction_type = Array.concat(junction_type,jtype);
				junction_defects = Array.concat(junction_defects,defect_coder(0,phase,jtype,x,y,"enc"));
				changeValues(10,10,100);
				
				//JUNCTION ORIENTATION
					//REALLY SIMPLE... needs improvement
					getSelectionBounds(xo, yo, wS, hS);
					for(yJ = yo-1; yJ < yo+hS+1; yJ++){
						for(xJ = xo-1; xJ < xo+wS+1; xJ++){
							if(getPixel(xJ,yJ)==13){xj_orient = xJ; yj_orient = yJ;}
						}
					}
					xy_j = xy_coder(0,0,xy,"dec");
					j_ndp = NDP(0,-1,xj_orient-xy_j[0],yj_orient-xy_j[1]);
					junction_orientations = Array.concat(junction_orientations, j_ndp);
					
					
			}
		}
	}
	run("Select None"); changeValues(13,13,12); // v029... reverts pixels.
	print("Junctions:");
	Array.print(junctions);
		
	// Collect Defects
	defects_array = Array.concat(defects_array,terminal_defects);
	defects_array = Array.concat(defects_array,junction_defects);
	print("TD:" + terminal_defects.length);
	print("JD:" + junction_defects.length);
	print("AD:" + defects_array.length);
	// Collect Orientations
	defects_orientations = Array.concat(defects_orientations,terminal_orientations);
	defects_orientations = Array.concat(defects_orientations,junction_orientations);
	print("TO:" + terminal_orientations.length);
	print("JO:" + junction_orientations.length);
	print("AO:" + defects_orientations.length);
	
	// Positive
	if(i_phase==0){
		pos_t_defects = terminal_defects; outputTD("pos_t_defects.L",pos_t_defects.length);
		pos_j_defects = junction_defects; outputTD("pos_j_defects.L",pos_j_defects.length);
	}
	// Negative
	else {
		neg_t_defects = terminal_defects; outputTD("neg_t_defects.L",neg_t_defects.length);
		neg_j_defects = junction_defects; outputTD("neg_j_defects.L",neg_j_defects.length);
	}
		
		
	} //}}}
	
	// END OF SKELETON ANALYSES // }}}

	
	// CORRELATION FUNCTION & ORDER PARAMETER CALCS //{{{
	
	/// Condition for correlation
	selectImage(image_positive_skeleton); pos_skeleton_pixel_count = 0;
	for(y=0; y<h; y++){for(x=0; x<w; x++){ if(getPixel(x,y)==255){ pos_skeleton_pixel_count++;} } } 
	pos_coverage_metric = pos_skeleton_pixel_count * wfft_period_pixel / (w * h);
	outputTD("Skel.Coverage.Metric.pos",pos_coverage_metric);
	selectImage(image_negative_skeleton); neg_skeleton_pixel_count = 0;
	for(y=0; y<h; y++){for(x=0; x<w; x++){ if(getPixel(x,y)==255){ neg_skeleton_pixel_count++;} } } 
	neg_coverage_metric = neg_skeleton_pixel_count * wfft_period_pixel / (w * h);
	outputTD("Skel.Coverage.Metric.neg",neg_coverage_metric);
	if(pos_coverage_metric > 0.4  || neg_coverage_metric > 0.4){
	
	/**
	Changes relative to: 2014_OP_Angles_29_20140320_Loop_Tidy.ijm
	replace "save_subfolder" with "save_subfolder"; comment-out lines for making folder.
	replace "save_title" with "opa_save_title"
	redefine "opa_save_title"
	opa_set_ds = true; {downsampling on... this is absolutely necessary!}
	ds_random = true; {random is selected; better distribution than Fibonacci downsampling.
	ds_fibonacci = true;
	opa_smooth_fs = round(wfft_period_pixel/2); // symmetric median filter radius
	opa_smooth_ws = round(wfft_period_pixel/2); // window filter radius
	replace "followTwoEncode(" with "followTwoEncodeVS("
	include "_CL_" in save instructions
	save image
	close plots
	
	Functions Required (Crossover):
	o	xy_coder()
	o	setForegroundIndex()
	o	neighbourValueExact()
	o>	followTwoEncode() ** followTwoEncodeVS()
	>	components()
	>	random_selector()
	>	downsample()
	>	downsampleFibonacci()
	>	arraySmoothFilterSymmetric()
	>	arraySmoothWindowSymmetric()
	>	base62converter()
	Functions Unused (Unique?):
	>	arraySmoothFilter()
	>	arraySmoothWindow()
	**/
	
	if(pos_coverage_metric > 0.4){ selectImage(image_positive_skeleton); outputTD("Correlation.Phase",1);}
	else { selectImage(image_negative_skeleton); outputTD("Correlation.Phase",0);}
	run("Duplicate...", "title=positive_skeleton_correlation");
	positive_correlation_skeleton = getImageID();
	opa_save_title = "correlation_data";
	
	// Settings
	
	/**
	smoothing/filtering should be set according to period
	
	**/
	
	d_sample = true;
	opa_factor = 50; // factor for downsampling
	opa_set_ds_points = 3500; // alternative downsampling, by approx number of points
	opa_set_ds = true; // if true, use opa_set_ds_points instead of opa_factor
	ds_random = true;
	ds_fibonacci = true; // Turns on only if TRUE and ds_random = FALSE
	opa_LT_max = 1; // number of times the loop should be repeated
	opa_smooth_fs = round(wfft_period_pixel/2); // symmetric median filter radius
	opa_smooth_ws = round(wfft_period_pixel/2); // window filter radius
	opa_plot_multi = true;
	
	outputTD("opa_factor",opa_factor); outputTD("opa_set_ds",opa_set_ds); outputTD("opa_set_ds_points",opa_set_ds_points);
	//save_subfolder = getDirectory("Choose save location.");
	//opa_save_title_old = substring(getTitle(),0,indexOf(getTitle(),"."));
	//getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	//month += 1; timestring = "" + year + "" + month + "" + dayOfMonth + hour + "" + minute + "" + second;
	//opa_save_title = opa_save_title_old + "+" + base62converter(parseInt(timestring),true);
	
	
	changeValues(1,255,255);
	
	opa_original = 255;
	opa_junction = 249; // needs to be smaller than the minimum "regular" value
	opa_terminal = 254;
	opa_single = 0;
	opa_regular = 253;
	opa_regular_T = 252;
	opa_regular_J = 251;
	opa_regular_L = 250;
	
	w = getWidth();
	h = getHeight();
	
	opa_terminal_array = newArray(0);
	for(y=0; y<h; y++){
		for(x=0; x<w; x++){
			if(getPixel(x,y)==255){
				opa_nv = neighbourValue(x,y,200);
				if(opa_nv>2){ setPixel(x,y,opa_junction); }
				if(opa_nv==1){ setPixel(x,y,opa_terminal); opa_terminal_array = Array.concat(opa_terminal_array,xy_coder(x,y,0,"enc"));}
				if(opa_nv==0){ setPixel(x,y,opa_single); }
				if(opa_nv==2){ setPixel(x,y,opa_regular); }
			}
		}
		showProgress(y,h);
	}
	
	// gets rid of all lines that have a terminal point
	opa_terminal_array_single = newArray(0);
	for(m=0; m<opa_terminal_array.length; m++){
		xy = xy_coder(0,0,opa_terminal_array[m],"dec");
		if(getPixel(xy[0],xy[1])==opa_terminal){ opa_terminal_array_single = Array.concat(opa_terminal_array_single,opa_terminal_array[m]); } // array with unique starts
		doWand(xy[0],xy[1],1,"8-connected");
		changeValues(opa_regular,opa_terminal,opa_regular_T); //run("Select 
	}
	run("Select None");
	
	// this gets rid of all lines that are connected to junctions only
	setForegroundIndex(opa_regular_J);
	opa_jterminal_array_single = newArray(0);
	for(y=0; y<h; y++){
		for(x=0; x<w; x++){
			if(getPixel(x,y)==opa_regular){
				opa_nv = neighbourValueExact(x,y,opa_junction);
				if(opa_nv==1){
					floodFill(x,y,"8-connected");
					opa_jterminal_array_single = Array.concat(opa_jterminal_array_single,xy_coder(x,y,0,"enc"));
				}
			}
		}
		showProgress(y,h);
	}
	
	
	setForegroundIndex(opa_regular_L);
	opa_loop_array = newArray(0);
	for(y=0; y<h; y++){
		for(x=0; x<w; x++){
			if(getPixel(x,y)==opa_regular){
				opa_loop_array = Array.concat(opa_loop_array,xy_coder(x,y,0,"enc"));
				floodFill(x,y,"8-connected");
			}
		}
		showProgress(y,h);
	}
	
	
	angle_average_length = 15;
	
	opa_nonLoop_array = Array.concat(opa_jterminal_array_single,opa_terminal_array_single);
	
	all_angles_array = newArray(0);
	all_x_array = newArray(0);
	all_y_array = newArray(0);
	all_x_comps = newArray(0);
	all_y_comps = newArray(0);
	all_starts_array = newArray(0);
	
	for(m=0; m<opa_nonLoop_array.length; m++){
		// Follow & create array of x's & y's
		xym = xy_coder(0,0,opa_nonLoop_array[m],"dec");
		pxvalue = getPixel(xym[0],xym[1]);
		encoded = followTwoEncodeVS(xym[0],xym[1],pxvalue,pxvalue);
		opa_e_length = encoded.length;
		x_array = newArray(encoded.length); y_array = newArray(encoded.length); angles = newArray(encoded.length);
		if(encoded.length > 2){ // This is to avoid NaN errors for single-pixels
			for(q=0; q<encoded.length; q++){
				xy = xy_coder(0,0,encoded[q],"dec");
				x_array[q] = xy[0]; y_array[q] = xy[1];
			}
			// Calculate Angles & Smooth
			x_component = components(x_array,angle_average_length);
			y_component = components(y_array,angle_average_length);
			/** This is to ensure angles from 0 to 180 **/
			for(i=0; i<y_component.length; i++){
				if(y_component[i]<0){
					y_component[i] = -1*y_component[i];
					x_component[i] = -1*x_component[i];
				}
				angles[i] = atan(y_component[i] / x_component[i]);
			}
			/** this creates a purely array-based computation of correlation function */
			all_angles_array = Array.concat(all_angles_array,angles);
			all_x_array = Array.concat(all_x_array,x_array);
			all_y_array = Array.concat(all_y_array,y_array);
			all_x_comps = Array.concat(all_x_array,x_component);
			all_y_comps = Array.concat(all_y_array,y_component);
			if(m==0){ all_starts_array = Array.concat(0,opa_e_length);}
			else{ all_starts_array = Array.concat(all_starts_array,opa_e_length);}
			//print(opa_e_length);
		}
		showProgress(m,opa_nonLoop_array.length);
	}
	
	/** Herman's Orientation Parameter **/
	opa_x_component_sum = 0; opa_x_normalized = newArray(all_angles_array.length);
	opa_y_component_sum = 0; opa_y_normalized = newArray(all_angles_array.length);
	
	for(m=0; m<all_angles_array.length; m++){
		opa_x_normalized[m] = cos(2*all_angles_array[m]);
		opa_y_normalized[m] = sin(2*all_angles_array[m]);
		opa_x_component_sum += opa_x_normalized[m];
		opa_y_component_sum += opa_y_normalized[m];
	}
	
	// error if the sum is exactly zero for both
	opa_length_sum = sqrt(pow(opa_x_component_sum,2) + pow(opa_y_component_sum,2));
	opa_x_component_normalized = opa_x_component_sum / opa_length_sum;
	opa_y_component_normalized = opa_y_component_sum / opa_length_sum;
	opa_angle_sum = atan(opa_y_component_normalized / opa_x_component_normalized)  / 2; // due to 180 to 360 to 180 conversion
	opa_angle_sum_deg = opa_angle_sum * 180 / PI; 
	opa_length_normalized = opa_length_sum / all_angles_array.length;
	
	opa_phi_average = 0;
	for(m=0; m<all_angles_array.length; m++){
		opa_phi_average += pow(cos(all_angles_array[m]-opa_angle_sum),2);
	}
	opa_phi_average = opa_phi_average / all_angles_array.length;
	opa_Hermans = (3*opa_phi_average - 1) / 2;
	opa_pre_downsample_length = all_angles_array.length;
	print("Herman's OP: " + opa_Hermans); outputTD("opa_Hermans",opa_Hermans);
	
	print("X component:"+ opa_x_component_sum); print("Y component:"+ opa_y_component_sum);
	print("X comp norm:"+ opa_x_component_normalized); print("Y comp norm:"+ opa_y_component_normalized);
	print("Average Angle: " + opa_angle_sum_deg);
	
	// End of Herman's Order Parameter Calculation
	
	/** CORRELATION CALCULATION **/
	print("array length unsampled: "+all_x_array.length ); outputTD("Array_Length_initial",all_x_array.length);
	
	// LOOP FOR TESTING
	opa_loop_testing = true;
	if(opa_loop_testing){
	opa_aaa_Orig = Array.copy(all_angles_array);
	opa_axa_Orig = Array.copy(all_x_array);
	opa_aya_Orig = Array.copy(all_y_array);
	opa_ds_counts = newArray(opa_LT_max);
	opa_ds_times = newArray(opa_LT_max);
	opa_ds_CL = newArray(opa_LT_max);
	}
	
	opa_LT_count = 0;
	//opa_LT_max = 25;
	while(opa_LT_count<opa_LT_max){ // START OF LOOP TEST
		if(opa_loop_testing){
			all_angles_array = Array.copy(opa_aaa_Orig);
			all_x_array = Array.copy(opa_axa_Orig);
			all_y_array = Array.copy(opa_aya_Orig);
		}
	/** Moved these here for v24 ... otherwise cumulative! **/
	opa_diagonal = floor(sqrt(w*w + h*h))+2; //2 because length doesn't include 0 & floor-ing.
	opa_counts = newArray(opa_diagonal); // how many times a measurement is for a given distance
	opa_sums = newArray(opa_diagonal); // total of the correlation.
		
	//d_sample = false;
	//opa_factor = 100;
	//ds_random = true;
	if(opa_set_ds){if(opa_pre_downsample_length<=opa_set_ds_points){d_sample = false;}}
	if(d_sample){
		if(opa_set_ds){opa_factor = round(all_x_array.length/opa_set_ds_points);}
		//print("Down-Sample: ON!");
		if(ds_random){
			ds_i = random_selector(opa_factor,all_x_array.length);
			opa_aaa = newArray(ds_i.length); opa_axa = newArray(ds_i.length); opa_aya = newArray(ds_i.length);
			for(m=0; m<ds_i.length; m++){
				opa_aaa[m] = all_angles_array[ds_i[m]];
				opa_axa[m] = all_x_array[ds_i[m]];
				opa_aya[m] = all_y_array[ds_i[m]];
			}
			all_angles_array = opa_aaa;
			all_x_array = opa_axa;
			all_y_array = opa_aya;
		}
		else if(ds_fibonacci){
			all_angles_array = downsampleFibonacci(all_angles_array,opa_factor);
			all_x_array = downsampleFibonacci(all_x_array,opa_factor);
			all_y_array = downsampleFibonacci(all_y_array,opa_factor);
		}
		else{
			all_angles_array = downsample(all_angles_array,opa_factor);
			all_x_array = downsample(all_x_array,opa_factor);
			all_y_array = downsample(all_y_array,opa_factor);
		}
	}
	print("array length downsampled: "+all_x_array.length); opa_post_downsample_length = all_x_array.length;
	outputTD("Array_Length_downsampled",opa_post_downsample_length);
	
	opa_period = 0; // this is introduced to avoid self-sampling & over-weighting of neighbour points
	m_items = all_angles_array.length;
	total = m_items * (m_items - opa_period - 1) / 2;
	completion = 0;
	showStatus("Completion: " + round(100*completion/total));
	opa_time_start = getTime();
	for(ma=0; ma<m_items-1; ma++){
		angle_a = all_angles_array[ma];
		xa = all_x_array[ma]; ya = all_y_array[ma];
		for(mb = ma+1+opa_period; mb<m_items; mb++){
			/** Alternat fix: could change the ma+1 offset to something larger; greater than period?**/
			angle_b = all_angles_array[mb];
			corr = cos( 2*(angle_a - angle_b) );
			dist = round( sqrt( pow(xa-all_x_array[mb],2) + pow(ya-all_y_array[mb],2) ) );
			opa_sums[dist] += corr; opa_counts[dist] += 1;
		}
		completion += m_items - ma - 1;
		showStatus("Completion: " + round(100*completion/total));
		showProgress(ma,m_items);
	}
	
	opa_time_end = getTime();
	opa_time_total = opa_time_end - opa_time_start;
	//print("Time: " + round(opa_time_total/1000) + " s");
	
	// Correlation Results
	opa_Correlation_temp = newArray(opa_sums.length); opa_Distance_axis = newArray(opa_sums.length);
	for(m=0; m<opa_sums.length; m++){
		if(opa_counts[m]==0){ opa_Correlation_temp[m] = 0; }
		else { opa_Correlation_temp[m] = opa_sums[m] / opa_counts[m]; }
		setResult("Sum",m,opa_sums[m]);
		setResult("Counts",m,opa_counts[m]);
		opa_Distance_axis[m] = m;
	}
	opa_Correlation_temp[0] = 1; // at distance = 0, the correlaiton is 1.
	opa_Correlation_smoothed_a = arraySmoothFilterSymmetric(opa_Correlation_temp,1,opa_smooth_fs);
	opa_Correlation_smoothed_b = arraySmoothWindowSymmetric(opa_Correlation_smoothed_a,1,opa_smooth_ws);
	
	one_over_e = exp(-1);
	find_CL = true; m = 0; mprev = 0;
	while(find_CL){
		value = opa_Correlation_smoothed_b[m];
		if(value<one_over_e){
			if(m==mprev){ correlation_length = 0; find_CL = false;}
			else{
				yi = value; yj = opa_Correlation_smoothed_b[mprev];
				xi = m; xj = mprev;
				slope = (yj-yi)/(xj-xi); intercept = yi - slope * xi;
				correlation_length = (one_over_e - intercept) / slope;
				find_CL = false;
			}
		}	
		else if(value==one_over_e){ correlation_length = m; find_CL = false;}
		mprev = m; m++;
	}
	
	// Correlation Length: Line-fit
	// Simple Linear Regression (without intercept)
	// only include correlation values > 0... b/c of log. 
	opa_CL_range = minOf(opa_Correlation_smoothed_b.length,round(correlation_length * 2)); // value of correlation fit strongly dependent on this. 
	opa_x_bar = 0; opa_xx_bar = 0; opa_xy_bar = 0; opa_yy_bar = 0; opa_bar_count = 0;
	for(m = 0; m<opa_CL_range; m++){
		yexp = opa_Correlation_smoothed_b[m];
		if(yexp>0){
			ylog = log(yexp); // note "log()" is actually natural log!
			opa_x_bar += m;
			opa_xx_bar = m*m;
			opa_xy_bar = m * ylog;
			opa_yy_bar = ylog * ylog;
			opa_bar_count++;
		}
	}
	opa_bar_slope = opa_xy_bar / opa_xx_bar;
	correlation_length_linear = -1 / opa_bar_slope;
	opa_bar_r_smaple_c_coeff = opa_xy_bar / sqrt( opa_xx_bar * opa_yy_bar );
	opa_bar_R_squared = opa_bar_r_smaple_c_coeff*opa_bar_r_smaple_c_coeff*1.0;
	// plot for fit:
	opa_Correlation_fit = newArray(opa_Correlation_smoothed_a.length);
	for(m=0; m<opa_Correlation_fit.length; m++){
		opa_Correlation_fit[m] = exp(-m/correlation_length_linear);
	}
	
	
	print("CORRELATION LENGTH (stepwise|line-fit|R2): "+ correlation_length +"," + correlation_length_linear + "," + opa_bar_R_squared);
	outputTD("correlation_length",correlation_length); outputTD("correlation_length_linear",correlation_length_linear); outputTD("opa_bar_R_squared",opa_bar_R_squared);
	//print("CORRELATION LENGTH (line-fit), px: "+ correlation_length_linear);
	//print("R-squared for line-fit: " + opa_bar_R_squared);
	
	if(opa_loop_testing){
	opa_ds_counts[opa_LT_count] = all_x_array.length;
	opa_ds_times[opa_LT_count] = opa_time_total;
	opa_ds_CL[opa_LT_count] = correlation_length;
	}
	
	// plots & results multi / opa_plot_multi
	if(opa_plot_multi){
		run("Clear Results");
		for(m=0; m<opa_sums.length; m++){
			setResult("Dist",m,opa_Distance_axis[m]);
			setResult("Correlation",m,opa_Correlation_temp[m]);
			setResult("C-Smooth-a",m,opa_Correlation_smoothed_a[m]);
			setResult("C-Smooth-b",m,opa_Correlation_smoothed_b[m]);
		}
		updateResults();
		selectWindow("Results"); saveAs("Measurements",save_subfolder+opa_save_title + "_CL_" + opa_LT_count + ".xls");
		
		Plot.create("Correlation Function " + opa_LT_count, "Distance (px)", "Correlation");
		Plot.setLimits(0, 400, -0.1, 1);
		Plot.setColor("red");
		Plot.add("circles", opa_Distance_axis,opa_Correlation_temp);
		Plot.setColor("magenta"); Plot.setLineWidth(3);
		Plot.add("dots", opa_Distance_axis,opa_Correlation_smoothed_a);
		Plot.setColor("blue");
		Plot.setLineWidth(2);
		Plot.add("line", opa_Distance_axis,opa_Correlation_smoothed_b);
		Plot.setColor("orange");
		Plot.setLineWidth(1);
		Plot.add("squares", opa_Distance_axis,opa_Correlation_fit);
		Plot.show();
		selectWindow("Correlation Function " + opa_LT_count); saveAs("PNG",save_subfolder+opa_save_title + "_CL_" + opa_LT_count + ".png"); close();
	} // end of opa_plot_multi for multi
	
	opa_LT_count++; // for LOOP TEST
	}//End of Correlation LOOP TESTING
	
	if(opa_loop_testing){
		print("\ndsCount,msTime,pxCL");
		for(m=0; m<opa_ds_CL.length; m++){
			print(opa_ds_counts[m]+","+ opa_ds_times[m]+","+opa_ds_CL[m]);
		}
	}
	
	if(opa_plot_multi){}else{
		run("Clear Results");
		for(m=0; m<opa_sums.length; m++){
			setResult("Dist",m,opa_Distance_axis[m]);
			setResult("Correlation",m,opa_Correlation_temp[m]);
			setResult("C-Smooth-a",m,opa_Correlation_smoothed_a[m]);
			setResult("C-Smooth-b",m,opa_Correlation_smoothed_b[m]);
		}
		updateResults();
		selectWindow("Results"); saveAs("Measurements",save_subfolder+opa_save_title + "_CL_" + opa_LT_count + ".xls");
		
		Plot.create("Correlation Function", "Distance (px)", "Correlation");
		Plot.setLimits(0, 400, -0.1, 1);
		Plot.setColor("red");
		Plot.add("circles", opa_Distance_axis,opa_Correlation_temp);
		Plot.setColor("magenta"); Plot.setLineWidth(3);
		Plot.add("dots", opa_Distance_axis,opa_Correlation_smoothed_a);
		Plot.setColor("blue");
		Plot.setLineWidth(2);
		Plot.add("line", opa_Distance_axis,opa_Correlation_smoothed_b);
		Plot.setColor("orange");
		Plot.setLineWidth(1);
		Plot.add("squares", opa_Distance_axis,opa_Correlation_fit);
		Plot.show();
		selectWindow("Correlation Function"); saveAs("PNG",save_subfolder+opa_save_title + "_CL_" + opa_LT_count + ".png"); close();
	} // end of opa_plot_multi for non-multi
	
	/** Save modified image **/
	selectImage(positive_correlation_skeleton); saveAs("PNG",save_subfolder+opa_save_title + "_skeleton_" + opa_LT_count + ".png"); close();
	
	//DOMAIN MAPPING //{{{
	//Source: 2014_OP_Angles_37_20140616_DomainMap.ijm
	//save_title = opa_save_title
	//domain_mapping = true;
	if(domain_mapping){
	
		run("Select None");

		selectImage(image_positive_lines);
		run("Duplicate...", "title=components_stack");
		components_stack = getImageID();
		changeValues(1,255,255);
		run("Make Binary");
		run("Distance Map");
		run("32-bit");
		pixel_ids = newArray(w*h);
		i = 0;
		for(y=0; y<h; y++){
			for(x=0; x<w; x++){
				pixel_ids[i] = getPixel(x,y);
				i++;
			}
		}
		Array.getStatistics(pixel_ids, pid_min, pid_max, pid_mean, pid_stdDev);
		changeValues(0,255,-10);
		
		
		a = 1/sqrt(2);
		b = 1;
		maxPx_d = -10;
		v1 = 0; v2 = 0; v3 = 0; v4 = 0; v5 = 0; v6 = 0; v7 = 0; v8 = 0;
		//setBatchMode(true);
		for(slice=0; slice<3; slice++){
			if(slice>0){run("Add Slice"); changeValues(-10,0,-10);}
			for(i=0; i<opa_aaa_Orig.length; i++){
				x = opa_axa_Orig[i];
				y = opa_aya_Orig[i];
				angle = opa_aaa_Orig[i];
				if(slice==0){setPixel(x,y,sin(angle/2));}
				if(slice==1){setPixel(x,y,cos(angle));}
				if(slice==2){setPixel(x,y,angle/abs(angle));}
			}
			
			pid = pid_max;
			while(pid>0){
				i = 0;
				for(y=0; y<h; y++){
					for(x=0; x<w; x++){
						if(pixel_ids[i]==pid && getPixel(x,y)==-10){
								       if(getPixel(x-1,y-1)>maxPx_d){n1=a; v1=getPixel(x-1,y-1);} else {n1=0;}
								       if(getPixel(x,y-1)>maxPx_d){n2=b; v2=getPixel(x,y-1);} else {n2=0;}
								       if(getPixel(x+1,y-1)>maxPx_d){n3=a; v3=getPixel(x+1,y-1);} else {n3=0;}
								       if(getPixel(x+1,y)>maxPx_d){	n4=b; v4=getPixel(x+1,y);} else {n4=0;}
								       if(getPixel(x+1,y+1)>maxPx_d){	n5=a; v5=getPixel(x+1,y+1);} else {n5=0;}
								       if(getPixel(x,y+1)>maxPx_d){	n6=b; v6=getPixel(x,y+1);} else {n6=0;}
								       if(getPixel(x-1,y+1)>maxPx_d){n7=a; v7=getPixel(x-1,y+1);	} else {n7=0;}
								       if(getPixel(x-1,y)>maxPx_d){n8=b; v8=getPixel(x-1,y);	} else {n8=0;}
								       V = v1*n1+v2*n2+v3*n3+v4*n4+v5*n5+v6*n6+v7*n7+v8*n8;
								       N = n1+n2+n3+n4+n5+n6+n7+n8;
								       if(N>0){
									 VN = V/N;
									 setPixel(x,y,VN);
								       }
						}
						i++;
					}
				}
				pid--;
			}
		showProgress(slice,3);
		} // end of slices
		
		
		w = getWidth();
		h = getHeight();
		pixel_ids = newArray(w*h);
		pid_phase = newArray(w*h);
		pid_x_comp = newArray(w*h);
		pid_y_comp = newArray(w*h);
		i = 0;
		for(y=0; y<h; y++){
			for(x=0; x<w; x++){
				if(getPixel(x,y)>-10){pixel_ids[i] = 1;}
				i++;
			}
		}
		i = 0;
		setSlice(1);
		for(y=0; y<h; y++){
			for(x=0; x<w; x++){
			if(pixel_ids[i]>0){pid_y_comp[i]=getPixel(x,y);}
			i++;
			}
		}
		
		i = 0;
		setSlice(2);
		for(y=0; y<h; y++){
			for(x=0; x<w; x++){
			if(pixel_ids[i]>0){pid_x_comp[i]=getPixel(x,y);}
			i++;
			}
		}
		
		i = 0;
		setSlice(3);
		for(y=0; y<h; y++){
			for(x=0; x<w; x++){
			if(pixel_ids[i]>0){pid_phase[i]=getPixel(x,y)/abs(getPixel(x,y));
			if(isNaN(pid_phase[i])){pid_phase[i]=1;}
			}
			i++;
			}
		}
		
		selectImage(image_positive_lines); 
		run("Duplicate...", "title=DomainMap");
		run("32-bit");
		//run("Add Slice");
		changeValues(-10,0,-10);
		i = 0;
		for(y=0; y<h; y++){
			for(x=0; x<w; x++){
			if(pixel_ids[i]>0){
			hyp = sqrt(pow(pid_x_comp[i],2)+pow(pid_y_comp[i],2));
			angle = acos(pid_phase[i]*pid_x_comp[i]/hyp);
			setPixel(x,y,angle);
			}
			i++;
			}
		}
		
		black = 0-PI/255;
		changeValues(-10,-10,black);
		LUT_angles(true);
		
		selectImage("DomainMap"); saveAs("TIFF",save_subfolder+"DomainMap.tif"); 
		//saveAs("TIFF",save_location+save_title + "_DomainMap" + ".tif");
		changeValues(255,255,0); // where the value should be zero, 255 pops up!
		saveAs("PNG",save_subfolder+"DomainMap.png"); close();
		// ADD ANGLE SCALE BA//{{{
			//Original: ADD_OPaS_Angle_Scale_Bar_01_20140516.ijm
		//}}} end of Add scalebar
		selectImage(components_stack); close();
		//setBatchMode(false);
		//print("DONE");
		
		
		function LUT_angles(true_or_false){
			if(true_or_false){
			reds = newArray(0, 254, 253, 252, 252, 251, 250, 249, 249, 248, 247, 246, 246, 245, 244, 244, 243, 242, 241, 241, 240, 239, 238, 238, 236, 234, 232, 230, 228, 226, 224, 222, 220, 218, 216, 214, 212, 209, 207, 205, 203, 201, 199, 197, 195, 193, 191, 187, 181, 175, 169, 163, 157, 151, 145, 139, 134, 128, 122, 116, 110, 104, 98, 92, 86, 80, 74, 69, 63, 57, 52, 50, 48, 45, 43, 41, 38, 36, 34, 31, 29, 27, 24, 22, 20, 18, 15, 13, 11, 8, 6, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 5, 7, 9, 11, 13, 15, 16, 18, 20, 22, 24, 26, 27, 29, 31, 33, 35, 37, 39, 40, 42, 45, 47, 50, 52, 54, 57, 59, 62, 64, 66, 69, 71, 74, 76, 79, 81, 83, 86, 88, 91, 93, 95, 98, 102, 106, 110, 114, 119, 123, 127, 131, 136, 140, 144, 149, 153, 157, 161, 166, 170, 174, 178, 183, 187, 191, 196, 199, 201, 204, 206, 208, 211, 213, 216, 218, 221, 223, 225, 228, 230, 233, 235, 237, 240, 242, 245, 247, 250, 252, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253);
			greens = newArray(0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 8, 9, 9, 11, 13, 14, 16, 18, 20, 22, 24, 26, 27, 29, 31, 33, 35, 37, 39, 40, 42, 44, 46, 48, 50, 51, 54, 56, 58, 61, 63, 66, 68, 71, 73, 75, 78, 80, 83, 85, 87, 90, 92, 95, 97, 99, 102, 104, 107, 109, 112, 115, 118, 121, 124, 127, 130, 132, 135, 138, 141, 144, 147, 150, 153, 156, 158, 161, 164, 167, 170, 173, 176, 176, 176, 176, 176, 177, 177, 177, 177, 177, 178, 178, 178, 178, 179, 179, 179, 179, 179, 180, 180, 180, 180, 180, 181, 181, 181, 181, 181, 181, 181, 181, 181, 181, 181, 182, 182, 182, 182, 182, 182, 182, 182, 182, 182, 182, 182, 184, 186, 188, 190, 192, 194, 196, 198, 199, 201, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221, 223, 225, 227, 227, 227, 226, 225, 225, 224, 223, 223, 222, 221, 221, 220, 219, 219, 218, 217, 217, 216, 216, 215, 214, 214, 213, 211, 207, 204, 200, 196, 192, 188, 185, 181, 177, 173, 169, 166, 162, 158, 154, 150, 147, 143, 139, 135, 131, 128, 124, 121, 118, 115, 112, 110, 107, 104, 101, 98, 95, 92, 89, 86, 84, 81, 78, 75, 72, 69, 53, 51, 31, 15);
			blues = newArray(0, 3, 6, 10, 13, 17, 20, 24, 27, 31, 34, 37, 41, 44, 48, 51, 55, 58, 62, 65, 69, 72, 75, 79, 82, 85, 88, 91, 94, 97, 100, 102, 105, 108, 111, 114, 117, 120, 123, 126, 129, 132, 135, 138, 141, 144, 146, 148, 150, 151, 153, 155, 156, 158, 159, 161, 162, 164, 165, 167, 168, 170, 171, 173, 174, 176, 177, 179, 180, 182, 183, 183, 183, 183, 183, 184, 184, 184, 184, 185, 185, 185, 185, 185, 186, 186, 186, 186, 186, 187, 187, 187, 187, 187, 187, 186, 185, 185, 184, 183, 182, 182, 181, 180, 180, 179, 178, 178, 177, 176, 176, 175, 174, 174, 173, 172, 171, 167, 163, 160, 156, 152, 148, 144, 140, 136, 132, 128, 125, 121, 117, 113, 109, 105, 101, 97, 94, 90, 86, 82, 78, 75, 71, 68, 64, 61, 57, 54, 50, 46, 43, 39, 36, 32, 29, 25, 22, 18, 15, 11, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
			setLut(reds,greens,blues);
			}
		}
		
		
		
		// NEED A MINIMUM LENGTH FUNCTION?
		
		/**
		Need a limit value for how many steps
		Should relate to the dimensions of the block copolymer
		3 x width of BCP?
		
		Could create a setting where if the LIMIT length is even, values get averaged differently.
		**/
	
		
	} //skip domain mapping
	//}}} end of DOMAIN MAPPING
	
	} // skip correlation
	
	// CORRELATION FUNCTION & ORDER PARAMETER CALCS //}}}

	
	// LINE-WIDTH ROUGHNESS or LINE-EDGE ROUGHNESS //{{{
	/** current version: LER-LWR_from_AF049c_integrating_05q_06.ijm **/
	/**  source: LER_new_multi_lines_06_20140124.ijm **/
	LER_analysis_option = 1;
	if(LER_analysis_option){

	// Settings //{{{

	min_area = positive_width*positive_width*10; // connected
	period = wfft_period_pixel; // pixels
	
	line_width = positive_width; //connected
	maxStretch = 2*line_width;
	oradius = line_width*1.5;
	
	initial_limit = wfft_period_pixel;
	
	// End of Settings //}}}
		
	// STEP 1: Remove Junctions & Edge-touching stretches
	/** 	Macro code to break all "lines"
		containing junctions into single,
		junctionless lines.
		Pre-requisite for measuring LER	**/
	
	// Usual Stuff
	selectImage("image_positive_skeleton"); // just for w & h purposes
	w = getWidth(); h = getHeight();
	
	// "optional"... so it doesn't need to be repeated.
	optional = 1;  //{{{
	if(optional){
	
	// Make "Broken" Images
	selectImage("image_positive_skeleton");
	image_positive_skeleton = getImageID();
	selectImage("image_positive_lines");
	image_positive_lines = getImageID();
	run("Duplicate...", "title=image_positive_broken");
	image_positive_broken = getImageID();
	
	
	// EDGE-TOUCHING STRETCHES
	// _edge_touching_stretches_3_20140123.ijm
	
	run("Select None"); run("Make Binary"); //just in case
	run("Analyze Particles...", "size=1-Infinity circularity=0.00-1.00 show=Nothing display clear record");
	
	// Semi-Circle Masks
	// Using a circular masks avoids most sharp corners
	// also less possibility of breaking lines... may need to check for that as a possibility...
	for(n=0; n<nResults; n++){
		x = getResult("XStart",n); y = getResult("YStart",n);	
		doWand(x,y,0,"8-connected");
		getSelectionCoordinates(xpoints,ypoints);
		xpoints = Array.concat(xpoints,xpoints[0]);
		ypoints = Array.concat(ypoints,ypoints[0]);
		x_touch = newArray(0); y_touch = newArray(0);
		for(i=0; i<xpoints.length-1; i++){
			check = false;
			if(xpoints[i]==0 && xpoints[i+1]==0){check = true; wbox = 2*oradius; hbox = abs(ypoints[i+1]-ypoints[i]); xbox = 0-oradius; ybox = ypoints[i+1];}  // x = 0
			if(ypoints[i]==0 && ypoints[i+1]==0){check = true; hbox = 2*oradius; wbox = abs(xpoints[i+1]-xpoints[i]); xbox = xpoints[i]; ybox = 0-oradius;}  // y = 0
			if(xpoints[i]==w && xpoints[i+1]==w){check = true; wbox = 2*oradius; hbox = abs(ypoints[i+1]-ypoints[i]); xbox = w -oradius; ybox = ypoints[i];}  // x = w
			if(ypoints[i]==h && ypoints[i+1]==h){check = true; hbox = 2*oradius; wbox = abs(xpoints[i+1]-xpoints[i]); xbox = xpoints[i+1]; ybox = h - oradius;}  // y = h
	
			
			if(check){
				xd = abs(xpoints[i]-xpoints[i+1]);
				yd = abs(ypoints[i]-ypoints[i+1]);
				length = maxOf(xd,yd);
				if(length>maxStretch){
					makeOval(xbox,ybox,wbox,hbox);
					changeValues(255,255,0);
				}
			
			}
		}
	}
	
	
	
	
	
	// Also need to remove any full loops
	
	// Break Lines
	w = getWidth(); h = getHeight();
	selectImage(image_positive_skeleton);
	setForegroundIndex(250);
	break_j_x = newArray(0);
	break_j_y = newArray(0);
	for(y=0; y<h; y++){
		for(x=0; x<w; x++){
			if(getPixel(x,y)>=250){
				nv = neighbourValue(x,y,250);
				if(nv>2){
					break_j_x = Array.concat(break_j_x,x);
					break_j_y = Array.concat(break_j_y,y);
				}
				if(nv!=2){ floodFill(x,y,"8-connected");}
			}
		}
	}
	// REMOVE LOOPS v05 fix
	loop_count = 0; setForegroundIndex(240);
	for(y=0; y<h; y++){
		for(x=0; x<w; x++){
			if(getPixel(x,y)==255){
				floodFill(x,y,"8-connected");
				loop_count++; //print("loop found");
				break_j_x = Array.concat(break_j_x,x);
				break_j_y = Array.concat(break_j_y,y);
			}
		}
	}
	changeValues(240,255,255); print("loop_count: " + loop_count); outputTD("loop_count",loop_count);
	
	
	selectImage(image_positive_broken);
	
	for(i=0; i<break_j_x.length; i++){
		x = break_j_x[i]; y = break_j_y[i];
		makeOval(x-period/2, y-period/2, period, period);
		changeValues(255,255,0);
	}
	run("Select None");
	
	// Analyze Particle Sizes
	// need to set measurements?
	run("Analyze Particles...", "size=1-Infinity circularity=0.00-1.00 show=Nothing display clear record");
	
	for(n=0; n<nResults; n++){
		area = getResult("Area",n);
		x = getResult("XStart",n); y = getResult("YStart",n);
		if(area<min_area){
			doWand(x,y,0,"8-connected");
			changeValues(255,255,0);
		}
		else{
			doWand(x,y,0,"4-connected"); // if anything is not 4-connected, remove it
			changeValues(255,255,254);
		}
	}
	run("Select None");
	changeValues(255,255,0);
	changeValues(254,254,255);
	
	// Slight Smoothing of Edges:
	for(y=0; y<h; y++){
		for(x=0; x<w; x++){
			if(getPixel(x,y)==255){
				nv = neighbourValue(x,y,254);
				if(nv<3){ setPixel(x,y,254); }
			}
		}
	}
	changeValues(254,254,0);
	
	//Broken Skeleton
	run("Select None");
	imageCalculator("Multiply create", "image_positive_broken","image_positive_skeleton");
	rename("image_broken_skeleton");
	image_broken_skeleton = getImageID();
	// Marking Skeleton
	for(y=0; y<h; y++){
		for(x=0; x<w; x++){
			if(getPixel(x,y)==255){
				nv = neighbourValue(x,y,1);
				if(nv==1){ setPixel(x,y,3); }
				if(nv==2){ setPixel(x,y,2); }
			}
		}
	}
	
	selectImage(image_positive_broken);
	run("Select None");
	// Create marked image
	/** Fix the problem by removing particles w/o lines inside * overly-small particles **/
	imageCalculator("Subtract create", "image_positive_broken","image_broken_skeleton");
	rename("image_broken_marked");
	image_broken_marked = getImageID();
	setThreshold(20, 255); run("Set Measurements...", "area mean standard min perimeter bounding shape feret's integrated redirect=None decimal=4");
	ipb_analysis = true;
	ipb_min_area = min_area; //Needs a better basis !@#$
	while(ipb_analysis){
		ipb_analysis = false;
		run("Analyze Particles...", "size=1-Infinity circularity=0.00-1.00 show=Nothing display clear record");
		for(n=0; n<nResults; n++){
			if(getResult("Mean",n)==255){ ipb_analysis = true;
				doWand(getResult("XStart",n),getResult("YStart",n),1,"8-connected");
				changeValues(1,255,10); run("Select None");
			}
			else if(getResult("Area",n)<ipb_min_area){ ipb_analysis = true;
				doWand(getResult("XStart",n),getResult("YStart",n),1,"8-connected");
				changeValues(1,255,15); run("Select None");
			}
		}
	}
	resetThreshold();
	/** End of Fix **/
	max_terminal_per = 0;
	for(n=0; n<nResults; n++){
		x = getResult("XStart",n); y = getResult("YStart",n);
		xy = findFirst(x,y,1,253);
		doWand(xy[0],xy[1],1,"8-connected");
		selectionCoordinatesExpanded();
		getSelectionCoordinates(xpoints,ypoints);
		//xskel = newArray(xpoints.length); Array.fill(xskel,-1);
		//yskel = newArray(xpoints.length); Array.fill(yskel,-1);
		count = 0;
		broken_terminals = newArray(0);
		for(i=0; i<xpoints.length; i++){
			x = xpoints[i]; y = ypoints[i]; value = getPixel(x,y);
			if(value==253 || value==252){
				setPixel(x,y,value-10); count++;
				if(value==252){bterm = xy_coder(x,y,0,"enc"); broken_terminals = Array.concat(broken_terminals,bterm);} //
				//if(CheckInArrayXY(xskel,yskel,x,y)){xskel[j] = x; yskel[j] = y;}
			}
		}
		if(broken_terminals.length==0){broken_terminals = newArray(2); broken_terminals[0] = xy_coder(xy[0],xy[1],0,"enc"); broken_terminals[1] = xy_coder(xy[0],xy[1],0,"enc"); } //v037
	
		max_terminal_per = maxOf(max_terminal_per,broken_terminals.length);
		setResult("Skel.Pts",n,count);
		setResult("T.Count",n, broken_terminals.length);
		/** if the break is done poorly, points may remain w/o terminals.. //!@#$
		May require a check-step where skeleton is subtracted; any lower than value are erased...**/
		if(broken_terminals.length<2){i_BT = broken_terminals.length;} else {i_BT = 2; }
		for(i=0; i<i_BT; i++){
			label = "Term"+i;
			setResult(i,n,broken_terminals[i]);
		}
	}
	run("Select None");
	//print(max_terminal_per);
	updateResults();
	
	
	changeValues(255,255,240); //required becaue followTwo uses pixel>=value
	
	}//}}} end of optional
	

	// MEASUREMENT OF SKELETON LENGTHS //{{{
	// This is, technically, unnecessary
	// But it does calculate the actual length of the skeleton!
	// So relevant to the final answer...
	for(n=0; n<nResults; n++){
		encoded = getResult("0",n);
		xy = xy_coder(0,0,encoded,"dec");
		dist = followTwo(xy[0],xy[1],241);
		setResult("Skel.Dist",n,dist[2]); updateResults(); //print("n: " + n);
		showStatus("Measuring Skeleton Lengths");
		showProgress(n,nResults);
	} //}}}
	updateResults();

	// ACTUAL MEASUREMENT OF LINE-EDGE DISPLACEMENTS  //{{{
	all_displacement_points_edge_x = newArray(0);
	all_displacement_points_edge_y = newArray(0);
	all_displacement_points_skel_x = newArray(0);
	all_displacement_points_skel_y = newArray(0);
	all_displacement_points_slope = newArray(0);
	all_displacement_points_width = newArray(0);
	
	
	displacement_points_edge_x = newArray(0);
	displacement_points_edge_y = newArray(0);
	displacement_points_skel_x = newArray(0);
	displacement_points_skel_y = newArray(0);
	displacement_points_slope = newArray(0);
	displacement_points_width = newArray(0);
	
	// LOOP For every PARTICLE //{{{
	for(np=0; np<nResults; np++){
		showStatus("Measuring LER & LWR");
		skel_time_i = getTime();
		
		// STEP 1: get the Skeleton, which is the centre line
		/** Should update to the new encoder **/
		encoded = getResult("0",np); // This is the first terminal point
		xy = xy_coder(0,0,encoded,"dec"); // decodes it	
		skeleton_points = followTwoEncode(xy[0],xy[1],241); // creates an encoded array
		xskel = newArray(skeleton_points.length);
		yskel = newArray(skeleton_points.length);
		for(i=0; i<skeleton_points.length; i++){
			xy = xy_coder(0,0,skeleton_points[i],"dec"); // decodes points
			xskel[i] = xy[0] + 0.5; // 0.5 places it in the actual centre.
			yskel[i] = xy[1] + 0.5; // 0.5 places it in the actual centre.
		}
		
		// Distance along skeleton:
		skel_distance = followTwoDistances(xy[0],xy[1],241);
		
		
		// Calculate Slope
		slope_skel = slopeArray(xskel,yskel,3,50);
		/** In future versions, try to smooth the skeleton points. **/
		
		skel_time_f = getTime();
		select_time_i = getTime();
		
		// STEP 2: get the Edge Points / Outline
			doWand(xy[0],xy[1],5,"8-connected"); // creates selection; requires some tolerance
			selectionCoordinatesExpanded(); // converts selection to pixel-by-pixel outer-edge basis
			getSelectionCoordinates(xpoints,ypoints); // records selection
			xpoints = Array.concat(xpoints,xpoints[0]); // to ensure a complete loop; no error from missing pair of points
			ypoints = Array.concat(ypoints,ypoints[0]); // to ensure a complete loop; no error from missing pair of points
				//polyEdgeSelectionInterpolate(xpoints,ypoints); // Shifts points to give a smoother edge
			//Reparameterizing from corners to mid-edge points to smooth the edges
			xtemp = newArray(xpoints.length);
			ytemp = newArray(xpoints.length);
			for(i=0; i<xpoints.length-1; i++){
				xtemp[i] = 0.5 * (xpoints[i] + xpoints[i+1]);
				ytemp[i] = 0.5 * (ypoints[i] + ypoints[i+1]);
			}
			xtemp[xtemp.length-1] = 0.5 * (xpoints[xtemp.length-1] + xpoints[0]);
			ytemp[ytemp.length-1] = 0.5 * (ypoints[ytemp.length-1] + ypoints[0]);
			xpoints = Array.copy(xtemp); ypoints = Array.copy(ytemp);
		//run("Select None");	//place this later... for visual effect
		setResult("Edge.Pts",np,xpoints.length); // how many edge points are there... give a sense of rejection rate!


		select_time_f = getTime();

		
		// Calculate the ROUGHNESS
		edge_time_i = getTime(); // figure out time-length correlation
		
		
		// Data Arrays for EDGE POINTS
		edge_points = xpoints.length;
		edge_array_valid = newArray(edge_points);
		
		edge_to_edge_distance = newArray(xpoints.length);
		edge_to_skel_distance = newArray(xpoints.length);
		edge_skel_avg_index = newArray(xpoints.length);
		skel_f_distance = newArray(xpoints.length);
		edge_array_cp = newArray(xpoints.length);		// Cross Product
		edge_array_closer = newArray(xpoints.length);
		edge_array_further = newArray(xpoints.length);
	
		// LOOP Over all EDGE Points //{{{
		range = 7;
		m_start = 0; m_end = xskel.length; m_range = range;
		i_start = 0; i_end = xpoints.length-1; i_range = range; //#!# need to give it a specific number
		for(n=0; n<xpoints.length; n++){ // loop over edge points; these are the start points
			xo = xpoints[n]; yo = ypoints[n];
			
			// EDGE -> SKELETON :1: find nearest points
			xc = -1; yc = -1; min_dist = 2*initial_limit;
			xm_dist_sq = 3*initial_limit; ym_dist_sq = 3*initial_limit;
			rerun  = true; 
			for(m=m_start; m<m_end; m++){ // loop over skeleton points to find the nearest point
				xi = xskel[m]; yi = yskel[m];
				x_dist_sq = (xo-xi)*(xo-xi); y_dist_sq = (yo-yi)*(yo-yi);
				if(x_dist_sq<xm_dist_sq || y_dist_sq < ym_dist_sq){
					dist = sqrt(x_dist_sq + y_dist_sq);
					//dist = sqrt( (xo-xi)*(xo-xi) + (yo-yi)*(yo-yi) );
					if(dist<min_dist){
						min_dist = dist; mm = m; //xc = xi; yc = yi; 
						xm_dist_sq = x_dist_sq; ym_dist_sq = y_dist_sq;
					}
				}
			}
			
			rerun = MinimaNotFoundRerunNecessary(m_start,m_end,mm,0,xskel.length);
			m_start = maxOf(0,mm-m_range); m_end = minOf(mm+m_range,xskel.length-1);
			if(rerun){ //skel_rerun_count++;
				for(m=0; m<xskel.length; m++){ // loop over skeleton points to find the nearest point
					xi = xskel[m]; yi = yskel[m];
					dist = sqrt( (xo-xi)*(xo-xi) + (yo-yi)*(yo-yi) );
					if(dist<min_dist){
						min_dist = dist; xc = xi; yc = yi; mm = m;
					}
				}
				m_start = maxOf(0,mm-m_range); m_end = minOf(mm+m_range,xskel.length-1);
				//if(mm_prev == mm){skel_pointless_rerun++;} //#!# for testing purposes
			}
			
			// EDGE -> SKELETON :2: Find the nearest point & angle at 90 degrees w.r.t. slope of skeleton:
			
			mm_start = maxOf(mm-2,0); mm_end = minOf(xskel.length-1,mm+2); ndp_pos = 2; ndp_neg = -2; mm_pos = -1; mm_neg = -1;
			for(m=mm_start; m<=mm_end; m++){
				ndp_m = NDP(1,slope_skel[m],xo-xskel[m],yo-yskel[m]);
				if(ndp_m>=0 && ndp_m<ndp_pos){ndp_pos = ndp_m; mm_pos = m;}
				else if(ndp_m<0 && ndp_m>ndp_neg){ndp_neg = ndp_m; mm_neg = m;}
			}
			
			// EDGE -> SKELETON :2: Correct difficult points & errors:
			if(mm_pos==-1){
				if(mm_neg==xskel.length-1){mm_pos = mm_neg-1;}
				else if(mm_neg==0){mm_pos = mm_neg+1; }
				else if(mm_neg==mm_start){mm_pos = maxOf(mm_neg-1,0);}
				else if(mm_neg==mm_end){mm_pos = minOf(mm_neg+1,xskel.length-1);}
				else {
					ndp_plus = NDP(1,slope_skel[mm_neg+1],xo-xskel[mm_neg+1],yo-yskel[mm_neg+1]);
					ndp_minus = NDP(1,slope_skel[mm_neg-1],xo-xskel[mm_neg-1],yo-yskel[mm_neg-1]);
					if(ndp_plus > ndp_minus){mm_pos = mm_neg+1;}
					else{mm_pos = mm_neg-1;}
				}
				ndp_pos = NDP(1,slope_skel[mm_pos],xo-xskel[mm_pos],yo-yskel[mm_pos]);
				// not sure if the assignment swap actually helps:
				if(ndp_pos < ndp_neg){ mm_temp = mm_neg; mm_neg = mm_pos; mm_pos = mm_temp; ndp_temp = ndp_neg; ndp_neg = ndp_pos; ndp_pos = ndp_temp;}
			}
			else if(mm_neg==-1){
				if(mm_pos==xskel.length-1){mm_neg = mm_pos-1;}
				else if(mm_pos==0){mm_neg = mm_pos+1;}
				else if(mm_pos==mm_start){mm_neg = maxOf(mm_pos-1,0);}
				else if(mm_pos==mm_end){mm_neg = minOf(mm_pos+1,xskel.length-1);}
				else {
					ndp_plus = NDP(1,slope_skel[mm_pos+1],xo-xskel[mm_pos+1],yo-yskel[mm_pos+1]);
					ndp_minus = NDP(1,slope_skel[mm_pos-1],xo-xskel[mm_pos-1],yo-yskel[mm_pos-1]);
					if(ndp_plus < ndp_minus){mm_neg = mm_pos+1; }
					else{mm_neg = mm_pos-1;}
				}
				ndp_neg = NDP(1,slope_skel[mm_neg],xo-xskel[mm_neg],yo-yskel[mm_neg]);
				// not sure if the assignment swap actually helps:
				if(ndp_pos < ndp_neg){ mm_temp = mm_neg; mm_neg = mm_pos; mm_pos = mm_temp; ndp_temp = ndp_neg; ndp_neg = ndp_pos; ndp_pos = ndp_temp;}
			}
	
			// Metrics:
			mm_pos_dist = sqrt(pow(xo-xskel[mm_pos],2)+pow(yo-yskel[mm_pos],2));
			mm_neg_dist = sqrt(pow(xo-xskel[mm_neg],2)+pow(yo-yskel[mm_neg],2));
	
			f_slope = (slope_skel[mm_pos]*mm_pos_dist + slope_skel[mm_neg]*mm_neg_dist)/(mm_pos_dist + mm_neg_dist);
			f = (xo + f_slope*yo - xskel[mm_neg] - f_slope*yskel[mm_neg])/(xskel[mm_pos]-xskel[mm_neg]+f_slope*(yskel[mm_pos]-yskel[mm_neg]));
			
			xc = f*xskel[mm_pos] + (1-f)*xskel[mm_neg]; // x position of point on skeleton / centre
			yc = f*yskel[mm_pos] + (1-f)*yskel[mm_neg]; // y position of point on skeleton / centre
			
			mm_metric = (f*mm_pos + (1-f)*mm_neg);
			// Unnecessary : All NDPs (normal dot products) will be 90 degrees.
				//ndp = NDP(1,f_slope,xpoints[n]-xs,ypoints[n]-ys);
				//ndp_degrees = acos(ndp)*180/PI;
				//print(round(ndp_degrees) + " "+ mm_metric);
	
			// EDGE -> SKELETON -> EDGE :: find intersection using first edge point, orthogonal skeleton point, and tranverse edge point
			if(mm_metric>=0 && mm_metric <= xskel.length-1){
				rerun = true;
				for(i=i_start; i<i_end; i++){ // re=loop over edge points to find the width
						vector = segmentIntersection(xo,yo,xc,yc,xpoints[i],ypoints[i],xpoints[i+1],ypoints[i+1]);
						if(vector[0]==1 && vector[2]<3.5){ // [2] constraint is to avoid hitting wrong side
							i_start = maxOf(0,i-i_range); i_end = minOf(i+i_range,xpoints.length-1);
							i = xpoints.length-1; rerun = false;
							xe = vector[4]; ye = vector[5];
						}
						else{xe = -1; ye = -1;}
				}
				if(rerun){ //edge_rerun_count++;
					for(i=0; i<xpoints.length-1; i++){ // re=loop over edge points to find the width
							vector = segmentIntersection(xo,yo,xc,yc,xpoints[i],ypoints[i],xpoints[i+1],ypoints[i+1]);
							if(vector[0]==1 && vector[2]<4){
								i_start = maxOf(0,i-i_range); i_end = minOf(i+i_range,xpoints.length-1);
								i = xpoints.length-1;
								xe = vector[4]; ye = vector[5];
							}
							else{xe = -1; ye = -1;}
					}
				}
	
			}
			else{xe = -1; ye = -1; vector = newArray(0,0,0,0,0,0);}
			// Record Data Points
			/** Could potentially do a second centre-to-tranverse-edge distance measurement **/
			if(xe==-1){ edge_to_edge_distance[n] = -1; edge_to_skel_distance[n] = -1; edge_array_cp[n] = 0;}
			else{	edge_array_valid[n] = 1;
				edge_to_edge_distance[n] = sqrt((xo-xe)*(xo-xe) + (yo-ye)*(yo-ye));
				edge_to_skel_distance[n] = sqrt((xo-xc)*(xo-xc) + (yo-yc)*(yo-yc));
				edge_array_cp[n] = nCrossProduv(xo-xc,yo-yc,1,f_slope);
			}
			skel_f_distance[n] = f*skel_distance[mm_pos] + (1-f)*skel_distance[mm_neg]; //print(f + " " + skel_distance[mm_neg] + " " + skel_distance[mm_pos]);
			edge_skel_avg_index[n] = (mm_pos + mm_neg)/2;
		} //}}} LOOP over all EDGE Points
		
		
		edge_time_f = getTime(); // figure out time-length correlation
		
		// CALCULATIONS & RESULTS //{{{
		// CALCULATE AVERAGE WIDTH (& STANDARD DEVIATIONS) FOR EACH LINE
		// FINALLY... the averages and the LER measurements
		calc_time_i = getTime();
		
		// POSITIVE (CP), NEGATIVE (CP), & Just VALID
		// LWR : WIDTHS / EDGE-to-EDGE distances
		p_edge_to_edge_distance = newArray(edge_points); p_eed = 0;
		n_edge_to_edge_distance = newArray(edge_points); n_eed = 0;
		v_edge_to_edge_distance = newArray(edge_points); v_eed = 0;
		// LER : EDGE-to-CENTRE/SKELETON distances
		p_edge_to_skel_distance = newArray(edge_points); p_esd = 0;
		n_edge_to_skel_distance = newArray(edge_points); n_esd = 0;
		v_edge_to_skel_distance = newArray(edge_points); v_esd = 0;
		
		for(m=0; m<edge_points; m++){
			if(edge_array_valid[m]==1 && edge_to_edge_distance[m]>0 && edge_to_skel_distance[m]>0){
				if(edge_to_edge_distance[m]<initial_limit*1.5 && edge_to_skel_distance[m]<initial_limit){ //points too far apart otherwise.
					if(edge_array_further[m]>1){ // old: edge_array_cp[m]==1
						p_edge_to_edge_distance[p_eed] = edge_to_edge_distance[m]; p_eed++; 
						p_edge_to_skel_distance[p_esd] = edge_to_skel_distance[m]; p_esd++;
					}
					if(edge_array_further[m]<=1){ // old: edge_array_cp[m]==-1
						n_edge_to_edge_distance[n_eed] = edge_to_edge_distance[m]; n_eed++; 
						n_edge_to_skel_distance[n_esd] = edge_to_skel_distance[m]; n_esd++;
					}
					v_edge_to_edge_distance[v_eed] = edge_to_edge_distance[m]; v_eed++; 
					v_edge_to_skel_distance[v_esd] = edge_to_skel_distance[m]; v_esd++;
				}
			}
			
		}
			//if(np==128){Array.print(v_edge_to_edge_distance);}
			//print(edge_points + " " + p_eed + " " + n_eed + " " + v_eed);
		
		// Trim Arrays
		p_edge_to_edge_distance = Array.trim(p_edge_to_edge_distance,p_eed);
		n_edge_to_edge_distance = Array.trim(n_edge_to_edge_distance,n_eed);
		v_edge_to_edge_distance = Array.trim(v_edge_to_edge_distance,v_eed);
		p_edge_to_skel_distance = Array.trim(p_edge_to_skel_distance,p_eed);
		n_edge_to_skel_distance = Array.trim(n_edge_to_skel_distance,n_eed);
		v_edge_to_skel_distance = Array.trim(v_edge_to_skel_distance,v_eed);
			//print(p_edge_to_edge_distance.length);
			//if(np==128){Array.print(v_edge_to_edge_distance);}
		// Get Statistics
		Array.getStatistics(p_edge_to_skel_distance, p_esd_min, p_esd_max, p_esd_mean, p_esd_stdev); // P Edge
		Array.getStatistics(n_edge_to_skel_distance, n_esd_min, n_esd_max, n_esd_mean, n_esd_stdev); // N Edge
		Array.getStatistics(v_edge_to_skel_distance, v_esd_min, v_esd_max, v_esd_mean, v_esd_stdev); // V Edge
		Array.getStatistics(p_edge_to_edge_distance, p_eed_min, p_eed_max, p_eed_mean, p_eed_stdev); // P Width
		Array.getStatistics(n_edge_to_edge_distance, n_eed_min, n_eed_max, n_eed_mean, n_eed_stdev); // N Width
		Array.getStatistics(v_edge_to_edge_distance, v_eed_min, v_eed_max, v_eed_mean, v_eed_stdev); // V Width
		
		calc_time_f = getTime();
		// Set Results //{{{
		setResult("P.E.sigma",np,p_esd_stdev); 
		setResult("P.E.avg",np,p_esd_mean); 
		setResult("P.E.count",np,p_esd); 
		setResult("P.E.sum",np,p_esd*p_esd_mean);
		setResult("P.E.min",np,p_esd_min);
		setResult("P.E.max",np,p_esd_max);
		setResult("P.W.sigma",np,p_eed_stdev); 
		setResult("P.W.avg",np,p_eed_mean); 
		setResult("P.W.count",np,p_eed); 
		setResult("P.W.sum",np,p_eed*p_eed_mean);
		setResult("P.W.min",np,p_eed_min);
		setResult("P.W.max",np,p_eed_max);
		
		setResult("N.E.sigma",np,n_esd_stdev); 
		setResult("N.E.avg",np,n_esd_mean); 
		setResult("N.E.count",np,n_esd); 
		setResult("N.E.sum",np,n_esd*n_esd_mean);
		setResult("N.E.min",np,n_esd_min);
		setResult("N.E.max",np,n_esd_max);
		setResult("N.W.sigma",np,n_eed_stdev); 
		setResult("N.W.avg",np,n_eed_mean); 
		setResult("N.W.count",np,n_eed); 
		setResult("N.W.sum",np,n_eed*n_eed_mean);
		setResult("N.W.min",np,n_eed_min);
		setResult("N.W.max",np,n_eed_max);
		
		setResult("V.E.sigma",np,v_esd_stdev); 
		setResult("V.E.avg",np,v_esd_mean); 
		setResult("V.E.count",np,v_esd); 
		setResult("V.E.sum",np,v_esd*v_esd_mean);
		setResult("V.E.min",np,v_esd_min);
		setResult("V.E.max",np,v_esd_max);
		setResult("V.W.sigma",np,v_eed_stdev); 
		setResult("V.W.avg",np,v_eed_mean); 
		setResult("V.W.count",np,v_eed); 
		setResult("V.W.sum",np,v_eed*v_eed_mean);
		setResult("V.W.min",np,v_eed_min);
		setResult("V.W.max",np,v_eed_max);
		
		
		setResult("t.skel",np,skel_time_f-skel_time_i);// figure out time-length correlation
		setResult("t.select",np,select_time_f-select_time_i);// figure out time-length correlation
		setResult("t.edge",np,edge_time_f-edge_time_i);// figure out time-length correlation
		setResult("t.calc",np,calc_time_f-calc_time_i);// figure out time-length correlation
		// End of Set Results //}}}
		// End of CALCULATIONS & RESULTS //}}}
		showProgress(n,nResults); run("Select None");
	} // End of loop//}}}
	updateResults(); //End of ACTUAL MEASUREMENT OF LINE-EDGE DISPLACEMENTS //}}}
	
	// CALCULATE AVERAGES //{{{
	// Averages are weighted by length
	P_length_total = 0; Psigma_length_total = 0;
	N_length_total = 0; Nsigma_length_total = 0;
	sigma_length_total = 0; length_total = 0;
	n_avg_length = 0; p_avg_length = 0; /** Define as zero... avoid problems: **/
	for(n=0; n<nResults; n++){
		Psigma = getResult("P.E.sigma",n);
		length = getResult("Skel.Dist",n);
		p_avg = getResult("P.E.avg",n); n_avg = getResult("N.E.avg",n);
		if( isNaN(Psigma) != 1 && Psigma > 0.1){  //Ensure that error values are not included
			Psigma_length_total += Psigma*length;
			P_length_total += length;
		}
		Nsigma = getResult("N.E.sigma",n);
		if( isNaN(Nsigma) != 1 && Nsigma > 0.1){  //Ensure that error values are not included
			Nsigma_length_total += Nsigma*length;
			N_length_total += length;
		}
		n_avg_length = n_avg*length; p_avg_length = p_avg*length; length_total += length;
	}
	Psigma_avg = Psigma_length_total / P_length_total; Nsigma_avg = Nsigma_length_total / N_length_total;
	print(Psigma_avg);
	
	outputTD("P.LER_sigma_avg",Psigma_avg);
	outputTD("N.LER_sigma_avg",Nsigma_avg);
	outputTD("LER_length_total_px",length_total);
	n_avg = n_avg_length / length_total; p_avg = p_avg_length / length_total; 
	outputTD("LER_N_width_avg",n_avg);
	outputTD("LER_P_width_avg",p_avg);
	outputTD("LER_width_avg",n_avg+p_avg);
	
	PNV_array = newArray("P","N","V");
	EW_array = newArray("E","W");
	sacsmm_array = newArray("sigma","avg","count","sum","min","max");
	
	for(pnv=0; pnv<PNV_array.length; pnv++){
		for(ew=0; ew<EW_array.length; ew++){
			for(sac=0; sac<sacsmm_array.length; sac++){
				items = newArray(nResults);
				column_label = PNV_array[pnv] + "." + EW_array[ew] + "." + sacsmm_array[sac];
				for(it=0; it<nResults; it++){
					items[it] = getResult(column_label,it);
				}
				Array.getStatistics(items,min,max,mean,sigma); count = nResults;
				outputTD(column_label+".avg",mean);
				outputTD(column_label+".sigma",sigma);
				outputTD(column_label+".min",min);
				outputTD(column_label+".max",max);
				outputTD(column_label+".count",count);
				print(column_label+".avg" + ": " + mean);
				print(column_label+".sigma" + ": " + sigma);
				print(column_label+".min" + ": " + min);
				print(column_label+".max" + ": " + max);
				print(column_label+".count" + ": " + count);
			}
		}
	}
	
	items = newArray(nResults);
	column_label = "Skel.Dist";
	for(it=0; it<nResults; it++){
		items[it] = getResult(column_label,it);
	}
	Array.getStatistics(items,min,max,mean,sigma); count = nResults;
	outputTD(column_label+".avg",mean);
	outputTD(column_label+".sigma",sigma);
	outputTD(column_label+".min",min);
	outputTD(column_label+".max",max);
	outputTD(column_label+".count",count);
	print(column_label+".avg" + ": " + mean);
	print(column_label+".sigma" + ": " + sigma);
	print(column_label+".min" + ": " + min);
	print(column_label+".max" + ": " + max);
	print(column_label+".count" + ": " + count);
	
	
	//End of CALCULATE AVERAGES //}}}
	
	// SAVE LER RESULTS //{{{
	selectWindow("Results");
	saveAs("Measurements",save_subfolder+"line_roughness_measurements.xls");
	run("Clear Results"); /**this needs to be turned back on**/
	//End of SAVE LER RESULTS //}}}
	
	} showStatus("LER & LWR Measurement Complete"); print("Done LER-LWR Test"); // END of LER Optional//}}}
	
	
	// RECORDING RESULTS //{{{
	
	for(i=0; i<defects_array.length; i++){
		values = defect_coder(defects_array[i],0,0,0,0,"dec");
		setResult("Phase",i,values[0]);
		setResult("Connectivity",i,values[1]);
		setResult("X",i,values[2]);
		setResult("Y",i,values[3]);
	}
	for(i=0; i<dot_defects.length; i++){
		values = defect_coder(dot_defects[i],0,0,0,0,"dec");
		setResult("Phase",i,values[0]);
		setResult("Connectivity",i,values[1]);
		setResult("X",i,values[2]);
		setResult("Y",i,values[3]);
	}
	updateResults();
	
	// END OF RECORDING RESULTS //}}}
	
	
	// CREATING COMPOSITES //{{{
	jn_radius = 0.5*minOf(positive_width,negative_width); edge_limit = maxOf(positive_width,negative_width); outputTD("DRAW_edge_limit",edge_limit); outputTD("DRAW_jn_radius",jn_radius);
	unusual_defects = true;
	regions_blocky = false;
	regions_fuzzy = false;
	regions_defects = true;
	regions_defects_new = true;
	figure_4 = true;
	
	// Defects + Regions //{{{
	if(regions_defects){
		selectImage(image004); // original-smoothed.
		run("Duplicate...", "title=yellow");
		changeValues(0,255,0);
		for(i=0; i<defects_array.length; i++){
			values = defect_coder(defects_array[i],0,0,0,0,"dec");
			phase = values[0]; connect = values[1]; x = values[2]; y = values[3];
			if(phase==0 && connect==1){ makeOval(x-jn_radius, y-jn_radius, jn_radius*2, jn_radius*2); changeValues(0,0,150);
				if(edgeDistance(w,h,x,y)<edge_limit){ changeValues(150,150,50); }}
			if(phase==1 && connect==1){ makeOval(x-jn_radius, y-jn_radius, jn_radius*2, jn_radius*2); changeValues(0,0,125);
				if(edgeDistance(w,h,x,y)<edge_limit){ changeValues(125,125,40); }}
			if(phase==0 && connect>=3){ makeOval(x-jn_radius, y-jn_radius, jn_radius*2, jn_radius*2); changeValues(0,0,250);}
			if(phase==1 && connect>=3){ makeOval(x-jn_radius, y-jn_radius, jn_radius*2, jn_radius*2); changeValues(0,0,200);}
		}
		run("Select None");
		selectImage(image003); // original-smoothed.
		run("Duplicate...", "title=grey");
		run("Multiply...", "value=0.5");
		selectImage(image_positive_lines);
		run("Duplicate...", "title=red");
		run("Multiply...", "value=0.5");
		selectImage(image_negative_lines);
		run("Duplicate...", "title=blue");
		run("Multiply...", "value=0.2");
		selectImage(image_positive_dots);
		run("Duplicate...", "title=magenta");
		run("Multiply...", "value=0.5");
		selectImage(image_negative_dots);
		run("Duplicate...", "title=cyan");
		run("Multiply...", "value=0.5");
		run("Merge Channels...", "c1=red c2=yellow c3=blue c4=grey c5=cyan c6=magenta create");
		rename("regions_defects");
		image_regions_defects = getImageID();
		selectImage(image_regions_defects); saveAs("TIF",save_subfolder+"image_regions_defects"+".tif"); saveAs("png",save_subfolder+"image_regions_defects"+".png"); //close();
	} //}}}
	
	// NEW Defects + Regions //{{{
	// Adapted to use branch symbols for junction-defects
	
	if(regions_defects_new){
		
		rdn_colours = newArray("red","green","blue");
		
		
		//selectImage(image003); // original-smoothed.
		//run("Duplicate...", "title=grey");
		//run("Multiply...", "value=0.14");
		
		//GREEN
		selectImage(image004); // original-smoothed.
		run("Duplicate...", "title=green");
		changeValues(0,255,0);
		
		//RED: Positive Lines
		selectImage(image_positive_lines);
		run("Duplicate...", "title=red");
		run("Multiply...", "value=1");
		
		// BLUE: Negative Lines
		selectImage(image_negative_lines);
		run("Duplicate...", "title=blue");
		run("Multiply...", "value=0.2");
		
		//MAGENTA: Positive Dots
		selectImage(image_positive_dots);  //**
		run("Duplicate...", "title=magenta");
		run("Multiply...", "value=1");
		
		//CYAN: Negative Dots
		selectImage(image_negative_dots); //**
		run("Duplicate...", "title=cyan");
		run("Multiply...", "value=0.7");
		
		
		
		// Marking Defects
		t1p = newArray(255,252,0); // T Positive
		j3p = newArray(255,209,98); // J-3 Positive
		j4p = newArray(255,180,5); // J-3 Positive
		j5p = newArray(255,160,130); // J-3 Positive
		
		t1n = newArray(0,235,255); // T Negative
		j3n = newArray(0,56,160); // J-3 Negative
		j4n = newArray(30,120,0); // J-3 Negative
		j5n = newArray(94,22,255); // J-3 Negative
		jn_scale = 1.5*jn_radius;
		
		for(i=0; i<defects_array.length; i++){
			angle = acos(defects_orientations[i]);
			values = defect_coder(defects_array[i],0,0,0,0,"dec");
			phase = values[0]; connect = values[1]; x = values[2]; y = values[3];
			if(phase==0 && connect==1){
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn_colours[rdn]);
					makeOval(x-jn_radius, y-jn_radius, jn_radius*2, jn_radius*2); changeValues(0,0,t1p[rdn]);
					//if(edgeDistance(w,h,x,y)<edge_limit){ changeValues(150,150,50); }
				}
			}
			if(phase==1 && connect==1){
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn_colours[rdn]);
					makeOval(x-jn_radius, y-jn_radius, jn_radius*2, jn_radius*2); changeValues(0,255,t1n[rdn]);
					//if(edgeDistance(w,h,x,y)<edge_limit){ changeValues(125,125,40); }
				}
			}
				
			// JUNCTIONS: Positive
			// J-3: TREFOIL
			if(phase==0 && connect==3){ 
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn_colours[rdn]);
					selectionCombo(Trefoil_x,Trefoil_y,angle,0,jn_scale,x,y,1); changeValues(0,255,j3p[rdn]); run("Select None");
				}
			}
			// J-4: PLUS
			if(phase==0 && connect==4){
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn_colours[rdn]);
					selectionCombo(Plus_x,Plus_y,angle,0,jn_scale,x,y,1); changeValues(0,255,j4p[rdn]); run("Select None");
				}
			}
			// J-5+: STAR
			if(phase==0 && connect>=5){
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn_colours[rdn]);
					selectionCombo(Star_x,Star_y,angle,0,jn_scale,x,y,1); changeValues(0,255,j5p[rdn]); run("Select None");
				}
			}
			
			// JUNCTIONS: Negative
			// J-3: TREFOIL
			if(phase==1 && connect==3){
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn_colours[rdn]);
					selectionCombo(Trefoil_x,Trefoil_y,angle,0,jn_scale,x,y,1); changeValues(0,255,j3n[rdn]); run("Select None");
				}
			}
			// J-4: PLUS
			if(phase==1 && connect==4){ 
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn_colours[rdn]);
					selectionCombo(Plus_x,Plus_y,angle,0,jn_scale,x,y,1); changeValues(0,255,j4n[rdn]); run("Select None");
				}
			}
			// J-5+: STAR
			if(phase==1 && connect>=5){
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn_colours[rdn]);
					selectionCombo(Star_x,Star_y,angle,0,jn_scale,x,y,1); changeValues(0,255,j5n[rdn]); run("Select None");
				}
			}
			
			
			
		}
		run("Select None");
		
		
		
		
		run("Merge Channels...", "c1=red c2=green c3=blue c5=cyan c6=magenta create"); // c4=grey
		rename("regions_defects_NEW");
		image_regions_defects_NEW = getImageID();
		selectImage(image_regions_defects_NEW); saveAs("TIF",save_subfolder+"image_regions_defects_NEW"+".tif"); saveAs("png",save_subfolder+"image_regions_defects_NEW"+".png"); //close();
	}//}}}
	
	
	// NEW Defects + Regions //{{{
	// Adapted to use branch symbols for junction-defects
	
	if(figure_4){
		
		rdn_colours = newArray("red","green","blue");
		rdn2_colours = newArray("redB","greenB","blueB");
		
		selectImage(image003); // original-smoothed.
		run("Duplicate...", "title=grey");
		run("Multiply...", "value=0.5");
		run("Duplicate...", "title=greyB");
		
		//GREEN
		selectImage(image004); // original-smoothed.
		run("Duplicate...", "title=green");
		changeValues(0,255,0);
		run("Duplicate...", "title=greenB");
		
		//RED: Positive Lines
		selectImage(image_positive_lines);
		run("Duplicate...", "title=red");
		run("Multiply...", "value=1");
		changeValues(0,255,0);
		run("Duplicate...", "title=redB");
		
		// BLUE: Negative Lines
		selectImage(image_negative_lines);
		run("Duplicate...", "title=blue");
		run("Multiply...", "value=0.2");
		changeValues(0,255,0);
		run("Duplicate...", "title=blueB");
		// The following is to prevent deletion of blueB from somewhere
		run("Duplicate...", "title=blueB2");
		
		//MAGENTA: Positive Dots
		selectImage(image_positive_dots);  //**
		run("Duplicate...", "title=magenta");
		run("Multiply...", "value=1");
		
		//CYAN: Negative Dots
		//selectImage(image_negative_dots); //**
		//run("Duplicate...", "title=cyanB");
		//run("Multiply...", "value=0.7");
		
		
		
		// Marking Defects
		t1p = newArray(255,252,0); // T Positive
		j3p = newArray(255,209,98); // J-3 Positive
		j4p = newArray(255,180,5); // J-3 Positive
		j5p = newArray(255,160,130); // J-3 Positive
		
		t1n = newArray(0,235,255); // T Negative
		j3n = newArray(0,56,160); // J-3 Negative
		j4n = newArray(30,120,0); // J-3 Negative
		j5n = newArray(94,22,255); // J-3 Negative
		jn_scale = 1.5*jn_radius;
		
		for(i=0; i<defects_array.length; i++){
			angle = acos(defects_orientations[i]);
			values = defect_coder(defects_array[i],0,0,0,0,"dec");
			phase = values[0]; connect = values[1]; x = values[2]; y = values[3];
			if(phase==0 && connect==1){
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn_colours[rdn]);
					makeOval(x-jn_radius, y-jn_radius, jn_radius*2, jn_radius*2); changeValues(0,0,t1p[rdn]);
					//if(edgeDistance(w,h,x,y)<edge_limit){ changeValues(150,150,50); }
				}
			}
			if(phase==1 && connect==1){
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn2_colours[rdn]);
					makeOval(x-jn_radius, y-jn_radius, jn_radius*2, jn_radius*2); changeValues(0,255,t1n[rdn]);
					//if(edgeDistance(w,h,x,y)<edge_limit){ changeValues(125,125,40); }
				}
			}
				
			// JUNCTIONS: Positive
			// J-3: TREFOIL
			if(phase==0 && connect==3){ 
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn_colours[rdn]);
					selectionCombo(Trefoil_x,Trefoil_y,angle,0,jn_scale,x,y,1); changeValues(0,255,j3p[rdn]); run("Select None");
				}
			}
			// J-4: PLUS
			if(phase==0 && connect==4){
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn_colours[rdn]);
					selectionCombo(Plus_x,Plus_y,angle,0,jn_scale,x,y,1); changeValues(0,255,j4p[rdn]); run("Select None");
				}
			}
			// J-5+: STAR
			if(phase==0 && connect>=5){
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn_colours[rdn]);
					selectionCombo(Star_x,Star_y,angle,0,jn_scale,x,y,1); changeValues(0,255,j5p[rdn]); run("Select None");
				}
			}
			
			// JUNCTIONS: Negative
			// J-3: TREFOIL
			if(phase==1 && connect==3){
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn2_colours[rdn]);
					selectionCombo(Trefoil_x,Trefoil_y,angle,0,jn_scale,x,y,1); changeValues(0,255,j3n[rdn]); run("Select None");
				}
			}
			// J-4: PLUS
			if(phase==1 && connect==4){ 
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn2_colours[rdn]);
					selectionCombo(Plus_x,Plus_y,angle,0,jn_scale,x,y,1); changeValues(0,255,j4n[rdn]); run("Select None");
				}
			}
			// J-5+: STAR
			if(phase==1 && connect>=5){
				for(rdn=0; rdn<3; rdn++){
					selectImage(rdn2_colours[rdn]);
					selectionCombo(Star_x,Star_y,angle,0,jn_scale,x,y,1); changeValues(0,255,j5n[rdn]); run("Select None");
				}
			}
			
			
			
		}
		run("Select None");
		
		
		
		
		run("Merge Channels...", "c1=red c2=green c3=blue c4=grey c6=magenta create"); // c4=grey
		rename("figure_4_b");
		image_figure_4_b = getImageID();
		selectImage(image_figure_4_b); saveAs("TIF",save_subfolder+"image_figure_4_b"+".tif"); saveAs("png",save_subfolder+"image_figure_4_b"+".png"); //close();
		
		//CYAN: Negative Dots
		selectImage(image_negative_dots); //**
		run("Duplicate...", "title=nayc");
		run("Multiply...", "value=0.7");
		
		run("Merge Channels...", "c1=redB c2=greenB c3=blueB c4=greyB c5=nayc create"); // c4=grey
		rename("figure_4_c");
		image_figure_4_c = getImageID();
		selectImage(image_figure_4_c); saveAs("TIF",save_subfolder+"image_figure_4_c"+".tif"); saveAs("png",save_subfolder+"image_figure_4_c"+".png"); //close();
	}//}}}
	
	
	// REGIONS COMPOSITES
	// Blocky //{{{
	if(regions_blocky){
		selectImage(image003); // original-smoothed.
		run("Duplicate...", "title=grey");
		run("Multiply...", "value=0.5");
		selectImage(image_positive_lines);
		run("Duplicate...", "title=red");
		run("Multiply...", "value=0.5");
		selectImage(image_negative_lines);
		run("Duplicate...", "title=blue");
		run("Multiply...", "value=0.2");
		selectImage(image_positive_dots);
		run("Duplicate...", "title=magenta");
		run("Multiply...", "value=0.5");
		selectImage(image_negative_dots);
		run("Duplicate...", "title=cyan");
		run("Multiply...", "value=0.5");
		run("Merge Channels...", "c1=red c3=blue c4=grey c5=cyan c6=magenta create");
		rename("regions_blocky");
		image_regions_blocky = getImageID();
		selectImage(image_regions_blocky); saveAs("png",save_subfolder+"image_regions_blocky"+".png"); //close();
	} //}}}
	// Fuzzy //{{{
	if(regions_fuzzy){
		selectImage(image003); // original-smoothed.
		run("Duplicate...", "title=grey");
		run("Multiply...", "value=0.5");
		selectImage(image_positive_lines);
		run("Duplicate...", "title=red");
		run("Multiply...", "value=0.5"); run("Gaussian Blur...", "sigma=2");
		selectImage(image_negative_lines);
		run("Duplicate...", "title=blue");
		run("Multiply...", "value=0.2"); run("Gaussian Blur...", "sigma=2");
		selectImage(image_positive_dots);
		run("Duplicate...", "title=magenta");
		run("Multiply...", "value=0.5"); run("Gaussian Blur...", "sigma=2");
		selectImage(image_negative_dots);
		run("Duplicate...", "title=cyan");
		run("Multiply...", "value=0.5"); run("Gaussian Blur...", "sigma=2");
		run("Merge Channels...", "c1=red c3=blue c4=grey c5=cyan c6=magenta create");
		rename("regions_fuzzy");
		image_regions_fuzzy = getImageID();
		selectImage(image_regions_fuzzy); saveAs("png",save_subfolder+"image_regions_fuzzy"+".png"); //close();
	} //}}}
	
	// VALIDATION IMAGES
	// Unusual Defects //{{{
	if(unusual_defects){
		jn_radius = 5;
		imageCalculator("Add create",positive_j_skeleton,negative_j_skeleton);
		unusual_temp = getImageID(); Ch2 = getTitle();
		changeValues(1,1,254); // Terminal Points
		changeValues(11,11,75); // 
		changeValues(12,12,50);
		for(i=0; i<defects_array.length; i++){
			values = defect_coder(defects_array[i],0,0,0,0,"dec");
			x = values[2]; y = values[3];
			if(values[1]<3){
				if(getPixel(x,y)!=254){ // this is contigent on above definitions
					//MARK
					makeOval(x-jn_radius, y-jn_radius, jn_radius*2, jn_radius*2);
					//shiftValues(0,0,80,1); // v029
					changeValues(0,0,80); // v030
					setPixel(x,y,255);
				}
			}
			else if(values[1]>3){
				//MARK
				makeOval(x-jn_radius, y-jn_radius, jn_radius*2, jn_radius*2);
				//shiftValues(0,0,120,1); // v029
				changeValues(0,0,120); // v030
				setPixel(x,y,255);
			}
		}
		run("Select None");
		selectImage(image_negative_skeleton); Ch3 = getTitle();
		selectImage(image_positive_skeleton); Ch1 = getTitle();
		run("Merge Channels...", "c1=&Ch1 c2=&Ch2 c3=&Ch3 create keep");
		image_validation_defects = getImageID();
		selectImage(image_validation_defects); saveAs("png",save_subfolder+"image_validation_defects"+".png"); //close();
		selectImage(unusual_temp); close();
	} //}}}
	
	// END OF CREATING COMPOSITES //}}}
	
	
	// FINAL CALCULATIONS & LOGGING //{{{
	// Categorize Defects
	edge_limit = 5; outputTD("edge_limit",edge_limit);
	
	p_dot_count = 0; p_edgedot_count = 0; p_edgeterminal_count = 0;
	p_terminal_count = 0; p_j3_count = 0;
	p_j4_count = 0; p_jx_count = 0; p_jx = 0;
	n_dot_count = 0; n_edgedot_count = 0; n_edgeterminal_count = 0;
	n_terminal_count = 0; n_j3_count = 0;
	n_j4_count = 0; n_jx_count = 0; n_jx = 0;
	for(i=0; i<defects_array.length; i++){
		values = defect_coder(defects_array[i],0,0,0,0,"dec");
		phase = values[0]; connect = values[1]; x = values[2]; y = values[3];
		if(phase==0){
			if(connect==1){if(edgeDistance(w,h,x,y)>edge_limit){ p_terminal_count++;}else{ p_edgeterminal_count++; }}
			if(connect==3){ p_j3_count++; }
			if(connect==4){ p_j4_count++; }
			if(connect>4){ p_jx_count++; p_jx += connect-2; }
			// dots
			//if(getResult("Dot",i)==1){if(getResult("OnEdge",i)==1){ p_edgedot_count++; } else { p_dot_count++; } }
		}
		if(phase==1){
			if(connect==1){if(edgeDistance(w,h,x,y)>edge_limit){ n_terminal_count++;}else{ n_edgeterminal_count++; }}
			if(connect==3){ n_j3_count++; }
			if(connect==4){ n_j4_count++; }
			if(connect>4){ n_jx_count++; n_jx += connect-2;}
			// dots
			//if(getResult("Dot",i)==1){if(getResult("OnEdge",i)==1){ n_edgedot_count++; } else { n_dot_count++; } }
		}
	
	}
	p_dot_count = pos_edge_dot_count; p_edgedot_count = pos_edge_dot_count; //little fix
	n_dot_count = neg_edge_dot_count; n_edgedot_count = neg_edge_dot_count; //little fix
	
	tab = "0000";
	print("\n"+"Value"+"\t"+"Pos"+"\t"+"Neg"+"\t"+"Tot"+"\t"+"Type");
	print("0"+"\t"+p_edgeterminal_count+"\t"+n_edgeterminal_count+"\t"+"Terminals-Edge"); outputTD("PTE",p_edgeterminal_count); outputTD("NTE",n_edgeterminal_count);
	print("1"+"\t"+p_terminal_count+"\t"+n_terminal_count+"\t"+"Terminals"); outputTD("PT",p_terminal_count); outputTD("NT",n_terminal_count);
	print("1"+"\t"+p_edgedot_count+"\t"+n_edgedot_count+"\t"+"Dots-Edge"); outputTD("PDE",p_edgedot_count); outputTD("NDE",n_edgedot_count);
	print("2"+"\t"+p_dot_count+"\t"+n_dot_count+"\t"+"Dots"); outputTD("PD",p_dot_count); outputTD("ND",n_dot_count);
	print("1"+"\t"+p_j3_count+"\t"+n_j3_count+"\t"+"Junctions"); outputTD("PJ3",p_j3_count); outputTD("NJ3",n_j3_count);
	print("2"+"\t"+p_j4_count+"\t"+n_j4_count+"\t"+"Junctions"); outputTD("PJ4",p_j4_count); outputTD("NJ4",n_j4_count);
	print("3+"+"\t"+p_jx_count+"\t"+n_jx_count+"\t"+"Junctions"); outputTD("PJx",p_jx_count); outputTD("NJx",n_jx_count);
	
	// Total Defects
	p_total_defects = 2*p_dot_count + 1*p_edgedot_count + 0*p_edgeterminal_count + 1*p_terminal_count + 1*p_j3_count + 2*p_j4_count + p_jx; outputTD("Ptot",p_total_defects);
	n_total_defects = 2*n_dot_count + 1*n_edgedot_count + 0*n_edgeterminal_count + 1*n_terminal_count + 1*n_j3_count + 2*n_j4_count + n_jx; outputTD("Ntot",n_total_defects);
	total_defects = p_total_defects + n_total_defects; outputTD("Total_Defects",total_defects);
	total_area_px = w*h; outputTD("Total_Area_px",total_area_px); 
	total_area_nm = total_area_px*nm_per_pixel*nm_per_pixel; outputTD("Total_Area_nm",total_area_nm); outputTD("Total_Area_um",total_area_nm/1000000);
	defect_pair_density_nm = total_defects / 2 / total_area_nm; outputTD("Defect_Density_nm",defect_pair_density_nm);
	defect_pair_density_um = defect_pair_density_nm * 1000000; outputTD("Defect_Density_um",defect_pair_density_um);
	
	print("\n"+"Total Pos. Defects: "+p_total_defects+"\n"+"Total Neg. Defects: "+n_total_defects+"\n"+"Total Defects: "+total_defects);
	print("\n"+"Total Area (px): "+total_area_px+"\n"+"Total Area (nm²): "+total_area_nm);
	print("\n"+"Defect Density (defects-pairs/nm²): "+defect_pair_density_nm+"\n"+"Defect Density (defects-pairs/µm²): "+defect_pair_density_um); 
	
	//END of F.C.&L.//}}}
	
	// Multi-Image Data Storage //{{{
	if(image_list.length>1){
		/** Insert Arrays Here */
		time_image_end = getTime() - time_zero;
		image_time = time_image_end - time_image_start;
		MULTI_defect_pair_density_um[img_i] = defect_pair_density_um;
		MULTI_time[img_i] = image_time;
	} //}}}
	
	// SAVE FILES //{{{
	selectWindow("Log");
	saveAs("Text", save_subfolder + image_log_name + ".txt");
	
	selectWindow("Results");
	saveAs("Measurements",save_subfolder+"defect_coordinates.xls");
	//run("Clear Results");
	
	// Save images:
	selectImage(image_positive_lines); title = getTitle(); saveAs("png", save_subfolder + title + ".png"); close();
	selectImage(image_positive_dots); title = getTitle(); saveAs("png", save_subfolder + title + ".png"); close();
	selectImage(image_negative_lines); title = getTitle(); saveAs("png", save_subfolder + title + ".png"); close();
	selectImage(image_negative_dots); title = getTitle(); saveAs("png", save_subfolder + title + ".png"); close();
	selectImage(image_positive_skeleton); title = getTitle(); saveAs("png", save_subfolder + title + ".png"); close();
	selectImage(image_negative_skeleton); title = getTitle(); saveAs("png", save_subfolder + title + ".png"); close();
	selectImage(positive_j_skeleton); title = getTitle(); saveAs("png", save_subfolder + title + ".png"); close();
	selectImage(negative_j_skeleton); title = getTitle(); saveAs("png", save_subfolder + title + ".png"); close();
	selectImage(image_positive_broken); title = getTitle(); saveAs("png", save_subfolder + title + ".png"); close();
	selectImage(image_broken_skeleton); title = getTitle(); saveAs("png", save_subfolder + title + ".png"); close();
	selectImage(image_broken_marked); title = getTitle(); saveAs("png", save_subfolder + title + ".png"); close();

	// DATA SUMMARY OUTPUT: outputTD + outputL
	run("Clear Results");
	/**  HORIZONTAL **/
	for(n=0; n<output_tags.length; n++){
		data = output_data[n];
		setResult(output_tags[n],0,data);
	}
	updateResults();
	selectWindow("Results");
	saveAs("Measurements",save_subfolder+"outputTD_horizontal.xls");
	
	run("Clear Results");
	/**  VERTICAL **/
	for(n=0; n<output_tags.length; n++){
		data = output_data[n];
		setResult("Label",n,output_tags[n]);
		setResult("Data",n,data);
	}
	updateResults();
	selectWindow("Results");
	saveAs("Measurements",save_subfolder+"outputTD_vertical.xls");

	// END OF SAVE FILES //}}}
	
	// CLOSE SUPERFLUOUS IMAGES //{{{
	if(image_list.length>1){ i = -1; n = 0; N = nImages;
		while(n<N){ if(isOpen(i)){ selectImage(i); close(); n++;} i--;}
	}//}}}

// NON-VIABILITY	
} else { //}}}
	// Non-Viable //{{{
	selectWindow("Log");
	print("\\Clear"); print("Not a viable image!");
	print(series_number + ": " + image_list[series_number]);
	saveAs("Text", save_folder + series_number + ".txt");
	} // Viable Image //}}}
	
	//!@#$ close & save the originally opened image.
	if(isOpen(image000)){selectImage(image000); close();}
} // ###### MULTI_IMAGE LOOP END ######  //}}}

	// Multi-Image Results Summary //{{{
	if(image_list.length>1){
		print("\\Clear");
		for(i=0; i<image_list.length; i++){
			print(i + " " + image_list[i] + "  " + MULTI_defect_pair_density_um[i] + "  " + round(MULTI_time[i]/1000));
		}
	}
//}}}

beep(); showStatus("done"); print("done"); //setBatchMode(false);
// END  //}}}





                                                                                               





