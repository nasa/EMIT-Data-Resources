# Getting Started with Earth Surface Mineral Dust Source Investigation (EMIT) data in Earthdata Search

---
**Note:** Users must have [NASA Earthdata Login](https://urs.earthdata.nasa.gov/home) to download data

## Step 1. Go to Earthdata Search and Login

Click **Earthdata Login** in the upper right portion of the Earthdata Search Landing Page

Enter your Earthdata Login credentials (username/password)  

![Earthdata Landing Page](https://i.imgur.com/CMzS6kA.jpeg)


## Step 2. Finding EMIT Data

1. In the top left poriton of the screen type in "EMIT" (quotes must be used in order to perform the search). This will search for all EMIT products the LP DAAC manages.

![Earthdata Landing Page](https://i.imgur.com/UeuG3bd.jpeg)

2. Currently the LP DAAC manages two EMIT collections.
EMIT L2A Estimated Surface Reflectance and Uncertainty and Masks 60m V001 and EMIT L1B At-Sensor Calibrated Radiance and Geolocation Data 60m V001.

![Earthdata Landing Page](https://i.imgur.com/j2zWRR2.png)

Select the collection of granules that you want to view. If you are unsure which collection to choose the ***View Collection Details*** offers more information about each data product. The current state of the collection includes a very limited data release with only a handful of scenes from select locations.

## Step 3. Download EMIT Data

To view the scenes and location, click on the image to focus it on the map the location is indicated by a green box. A preview of the scene is shown on the right.

1. If you are interested in downloading all EMIT scenes discovered by Earthdata search, select the ***Download All*** icon at the bottom of the search results page.

![Earthdata Landing Page](https://i.imgur.com/YjcAZvU.png)

2. If you are interested in only a singular scene, click the download arrow on the card (or list row) for the select scene. You will then be able to download the select files related to that image.

![Earthdata Landing Page](https://i.imgur.com/cdOsgm3.png)

You can also access S3 information (e.g., AWS region, bucket, temperorary credentials for S3 access, and file names)

![Earthdata Landing Page](https://i.imgur.com/Y4jVdD0.png)

3. To download multiple granules, click on the green + symbol to add files to the project. Click on the green button on the bottom of the page that says **Download**. This will take you to the options page for customizing the download.

![Earthdata Landing Page](https://i.imgur.com/ouESYWJ.png)

## Step 4. Data Access

On the next page click the **Direct Download** option and click the green **Download Data** on the bottom left side of the page.

![Earthdata Landing Page](https://i.imgur.com/NuoENO8.png)

## Step 5. Download Status

![Earthdata Landing Page](https://i.imgur.com/T3sbhau.png)

The final page of instructions show the download links for data access in the cloud. You should see four tabs: Download Files, AWS S3 Access, Download Script, and Browse Imagery.

The **Download Files** tab provides the *https://* links for downloading the files locally

The **AWS S3 Access** tab provides the  *s3://* links, which is what we would use to access the data directly in region (us-west-2) within the AWS cloud.


[def]: ScreenshotStep1.png