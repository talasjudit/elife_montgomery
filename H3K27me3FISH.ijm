run("Duplicate...", " ");
run("8-bit");
run("Gaussian Blur...", "sigma=3");
// Input a tissue threshold percentage 0-100
tissueThreshPerc = 97;
nBins = 256;
getHistogram(values, count, nBins);
size = count.length;
// find culmulative sum
cumSum = 0;
for (i = 0; i<count.length; i++)
{
  cumSum += count[i];
}
//totalNumberOfPixels = getWidth() * getHeight()
//print(totalNumberOfPixels)
tissueValue = cumSum * tissueThreshPerc / 100;
print(tissueValue);
// cumulative sum of before
cumSumValues = count;
for (i = 1; i<count.length; i++)
{
  cumSumValues[i] += cumSumValues[i-1];
}
// find tissueValue
for (i = 1; i<cumSumValues.length; i++)
{
  if (cumSumValues[i-1] <= tissueValue && tissueValue <= cumSumValues[i])
    {// output tissue threshold:
    print(i);
    setThreshold(i,255);}
}
setOption("BlackBackground", false);
run("Convert to Mask");
roiManager("Select", 0);
run("Set Measurements...", "area area_fraction redirect=None decimal=3");
roiManager("Measure");
roiManager("Delete");