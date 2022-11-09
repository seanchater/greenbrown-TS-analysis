

var aoi = control_farm5_area2;
var start_year = 2011;
var end_year = 2022;
print('start year', start_year);
print('end year', end_year);
var years = ee.List.sequence(start_year, end_year);
var index_name = "NDVI";
// define function to calculate a spectral index to segment with LT
var segIndex = function(img) {
    var index = img.normalizedDifference(['SR_B5', 'SR_B4'])                      // calculate normalized difference (NDVI) of band 4 and band 7 (B4-B3)/(B4+B3)
                   .select([0], [index_name])                                    // ...name the band
                   .set('system:time_start', img.get('system:time_start')); // ...set the output system:time_start metadata to the input image time_start otherwise it is null
    return index;
};


//cloud mask
function maskL8srClouds(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = image.select('QA_PIXEL');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}

function clip_image(image) {
  return image.clip(aoi);
}
function slc_correct(image) {
  var filled1a = image.focal_mean(3, 'square', 'pixels', 2)
  return filled1a.blend(image);
}
function get_year_image(year) {
  var start = year + "-01-01";
  var end = year+ "-12-12";
  var l7composite = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
                  .filterBounds(aoi)
                  .filterDate(start, end) //consider only filtering appropriate season images
                  .filterMetadata('CLOUD_COVER','less_than',25)
                  .select(['SR_B3', 'SR_B4'], ['red', 'nir'])
                  .map(clip_image)
                  .map(slc_correct)
                  .median();
                  // .map(maskL8srClouds);
  return l7composite;
}

// collect a median monthly composite from L8 that has been filtered for clouds
function get_month(year, month) {
  var start = year + '-'+ month+ "-01"; // start of the month
  var end = ee.Date(start).getRange('month').end(); // get end of the month by using date ee.Date().getRange('month')
  var l8composite = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
                  .filterBounds(aoi)
                  .filterDate(start, end) //consider only filtering appropriate season images
                  // .filterMetadata('CLOUD_COVER','less_than',25)
                  .select(['SR_B4', 'SR_B5', 'QA_PIXEL'], ['red', 'nir', 'QA_PIXEL'])
                  .map(clip_image)
                  .map(maskL8srClouds)
                  .max();
                  // .median();
  return l8composite;
}


//########################################################################################################
//##### ANNUAL SR TIME SERIES COLLECTION BUILDING FUNCTIONS ##### 
//########################################################################################################

//----- MAKE A DUMMY COLLECTOIN FOR FILLTING MISSING YEARS -----
var dummyCollection = ee.ImageCollection([ee.Image([0,0,0,0,0,0]).mask(ee.Image(0)).rename(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']).toFloat()]); // make an image collection from an image with 6 bands all set to 0 and then make them masked values


//------ L8 to L7 HARMONIZATION FUNCTION -----
// slope and intercept citation: Roy, D.P., Kovalskyy, V., Zhang, H.K., Vermote, E.F., Yan, L., Kumar, S.S, Egorov, A., 2016, Characterization of Landsat-7 to Landsat-8 reflective wavelength and normalized difference vegetation index continuity, Remote Sensing of Environment, 185, 57-70.(http://dx.doi.org/10.1016/j.rse.2015.12.024); Table 2 - reduced major axis (RMA) regression coefficients
var harmonizationRoy = function(oli) {
  var slopes = ee.Image.constant([0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949]);        // create an image of slopes per band for L8 TO L7 regression line - David Roy
  var itcp = ee.Image.constant([-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029]);     // create an image of y-intercepts per band for L8 TO L7 regression line - David Roy
  var y = oli.select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'],['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']) // select OLI bands 2-7 and rename them to match L7 band names
             .resample('bicubic')                                                          // ...resample the L8 bands using bicubic
             .subtract(itcp.multiply(10000)).divide(slopes)                                // ...multiply the y-intercept bands by 10000 to match the scale of the L7 bands then apply the line equation - subtract the intercept and divide by the slope
             .set('system:time_start', oli.get('system:time_start'));                      // ...set the output system:time_start metadata to the input image time_start otherwise it is null
  return y.toUint16();//toShort();                                                                       // return the image as short to match the type of the other data
};


//------ RETRIEVE A SENSOR SR COLLECTION FUNCTION -----
var getSRcollection = function(year, month, sensor, aoi) {
  // get a landsat collection for given year, day range, and sensor
  // var srCollection = ee.ImageCollection('LANDSAT/'+ sensor + '/C01/T1_SR') // get surface reflectance images
  //                     .filterBounds(aoi)                                  // ...filter them by intersection with AOI
  //                     .filterDate(year+'-'+startDay, year+'-'+endDay);    // ...filter them by year and day range
  
  var start = year + '-'+ month+ "-01"; // start of the month
  var end = ee.Date(start).getRange('month').end(); // get end of the month by using date ee.Date().getRange('month')
  var srCollection = ee.ImageCollection("LANDSAT/" + sensor + "/C02/T1_L2")
                  .filterBounds(aoi)
                  .filterDate(start, end) //consider only filtering appropriate season images
                  .map(clip_image)
                  .map(maskL8srClouds);
  
  
  
  // apply the harmonization function to LC08 (if LC08), subset bands, unmask, and resample           
  srCollection = srCollection.map(function(img) {
    var dat = ee.Image(
      ee.Algorithms.If(
        sensor == 'LC08',                                                  // condition - if image is OLI
        harmonizationRoy(img.unmask()),                                    // true - then apply the L8 TO L7 alignment function after unmasking pixels that were previosuly masked (why/when are pixels masked)
        img.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'])                   // false - else select out the reflectance bands from the non-OLI image
           .unmask()                                                       // ...unmask any previously masked pixels 
           .resample('bicubic')                                            // ...resample by bicubic 
           .set('system:time_start', img.get('system:time_start'))         // ...set the output system:time_start metadata to the input image time_start otherwise it is null
      )
    );
    
    // make a cloud, cloud shadow, and snow mask from fmask band
    var qa = img.select('QA_PIXEL');                                       // select out the fmask band
    var mask = qa.bitwiseAnd(8).eq(0).and(                                 // include shadow
               qa.bitwiseAnd(16).eq(0)).and(                               // include snow
               qa.bitwiseAnd(32).eq(0));                                   // include clouds
    
    // apply the mask to the image and return it
    return dat.mask(mask); //apply the mask - 0's in mask will be excluded from computation and set to opacity=0 in display
  });

  return srCollection; // return the prepared collection
};


//------ FUNCTION TO COMBINE LT05, LE07, & LC08 COLLECTIONS -----
var getCombinedSRcollection = function(year, month, aoi) {
    var lt5 = getSRcollection(year, month, 'LT05', aoi);       // get TM collection for a given year, date range, and area
    var le7 = getSRcollection(year, month, 'LE07', aoi);       // get ETM+ collection for a given year, date range, and area
    var lc8 = getSRcollection(year, month, 'LC08', aoi);       // get OLI collection for a given year, date range, and area
    var mergedCollection = ee.ImageCollection(lt5.merge(le7).merge(lc8)); // merge the individual sensor collections into one imageCollection object
    return mergedCollection;                                              // return the Imagecollection
};


//------ FUNCTION TO REDUCE COLLECTION TO SINGLE IMAGE PER YEAR BY MEDOID -----
/*
  LT expects only a single image per year in a time series, there are lost of ways to
  do best available pixel compositing - we have found that a mediod composite requires little logic
  is robust, and fast
  
  Medoids are representative objects of a data set or a cluster with a data set whose average 
  dissimilarity to all the objects in the cluster is minimal. Medoids are similar in concept to 
  means or centroids, but medoids are always members of the data set.
*/

// make a medoid composite with equal weight among indices
var maxMosaic = function(inCollection, dummyCollection) {
  
  // fill in missing years with the dummy collection
  var imageCount = inCollection.toList(1).length();                                                            // get the number of images 
  var finalCollection = ee.ImageCollection(ee.Algorithms.If(imageCount.gt(0), inCollection, dummyCollection)); // if the number of images in this year is 0, then use the dummy collection, otherwise use the SR collection
  
  // calculate median across images in collection per band
  var max = finalCollection.max();                                                                       // calculate the median of the annual image collection - returns a single 6 band image - the collection median per band

  return max;
};


//------ FUNCTION TO APPLY max COMPOSITING FUNCTION TO A COLLECTION -------------------------------------------
var buildMosaic = function(year, month, aoi, dummyCollection) {                                                                      // create a temp variable to hold the upcoming annual mosiac
  var collection = getCombinedSRcollection(year, month, aoi);  // get the SR collection
  var img = maxMosaic(collection, dummyCollection)                     // apply the medoidMosaic function to reduce the collection to single image per year by medoid 
              .set('system:time_start', (new Date(year,month,1)).valueOf());  // add the year to each medoid image - the data is hard-coded  1st day of the current year-month combo
  return ee.Image(img);                                                   // return as image object
};


//------ FUNCTION TO BUILD ANNUAL MOSAIC COLLECTION ------------------------------
var buildMosaicCollection = function(startYear, endYear, aoi, dummyCollection) {
  var imgs = [];                                                                    // create empty array to fill
  for (var year = startYear; year <= endYear; year++) {                                      // for each year from hard defined start to end build medoid composite and then add to empty img array
    for (var month = 1; month <= 12; month++) {
      var tmp = buildMosaic(year, month, aoi, dummyCollection);               // build the medoid mosaic for a given year
      imgs = imgs.concat(tmp.set('system:time_start', (new Date(year,month,1)).valueOf()));  // concatenate the annual image medoid to the collection (img) and set the date of the image - hard coded to the (year-month) that is being worked on for 1st of that monthq
    }
    
  }
  return ee.ImageCollection(imgs);                                                  // return the array img array as an image collection
};


// get image collection of all monthly maxs from combined L5 L7 L8 
var final_collection = buildMosaicCollection(start_year, end_year, aoi, dummyCollection);
print("initial final:", final_collection);
// *********************************
// *********************************


// create a list to hold all the  median images
// var img_list = [];
// print(img_list)
// for (var year = start_year; year <= end_year; year++) {
//   // get a single landsat image in each year with the median
//   // var image = get_year_image(year);
  
//   for (var month = 1; month <= 12; month++) {
//     // get a single landsat 8 img of a median/max etc. monthly
//     var image = get_month(year, month);    
//     img_list.push(image);
    
//   }

//   // print(image);
//   // add to the list of yearly images
//   // img_list.push(image);
// }

// print(img_list);


// add some imgs to the map to see if they work
// for (var i = 0; i < 10; i++) {
//   var img = img_list[i];
//   Map.addLayer(img, null, 'img'+i);
// }
// smoosh all the images in year_list into an ImageCollection
// var final_collection = ee.ImageCollection.fromImages(img_list);

// create a projection from band 3 for export
var projection = final_collection.first().select('SR_B4').projection().getInfo();

// apply the function to calculate the segmentation index 
var final_collection = final_collection.map(segIndex); // ...set the output system:time_start metadata to the input image time_start otherwise it is null


// select only the NDVI band for export
var final_ndvi = final_collection.map(function(img) {
  return img.select("NDVI");
});
print("final_ndvi:", final_ndvi);


// var final_banded = final_ndvi.toBands();
// print("banded:", final_banded);

// Export.image.toDrive({
//   image: final_banded,
//   description: 'ndvi_2005-2022',
//   crs: projection.crs,
//   crsTransform: projection.transform,
//   scale: 30,
//   region: aoi
// });

Map.addLayer(final_ndvi.first(), null, "first");
// below export function gotten from https://github.com/fitoprincipe/geetools-code-editor/wiki/Batch
var options = {
  "name": "ndvi{id}",
  "scale": 30,
  "region": aoi,
  "type": "float"
};
var batch = require('users/fitoprincipe/geetools:batch').Download.ImageCollection.toDrive;
batch(final_ndvi, "", options);

// to run all and confirm all tasks automatically: https://github.com/gee-hydro/gee_monkey

// Map.addLayer(aoi);
// Map.centerObject(aoi, 16);

