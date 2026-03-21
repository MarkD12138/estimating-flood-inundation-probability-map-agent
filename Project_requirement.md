**Estimating a flood inundation probability map**

Your assignment is to build a spatial predictive model in R that
estimates ‘the probability that an area will be inundated by a flood”. 
The flood you are modeling is based on one that took place in Calgary,
Alberta (Canada) in 2013
(https://www.calgary.ca/water/flooding/history-calgary.html). It was the
worst flood there in over 100 years, and the Director of Disaster
Planning in *CITY X* has asked you to try a machine learning approach to
estimate what would happen if a Calgary-type weather event took place
there.

Train and validate a spatial predictive model using data from Calgary.
Deploy model to predict a flood in a **CITY X**. Write a report to the
"Director of Disaster Planning" with your findings. You get to choose
which City X you'd like to use for your model (more on that below).

To be clear, to do this analysis, you will *copy the code* from the land
conservation markdown (predicting the probability of land conservation)
and create new independent variables for a model that predicts
probability of flood inundation. You can wrangle data in ArcGIS or R,
it’s your choice, but your modeling will take place in R.  You must
stick to "glm" logistic regression only. 

**This is a team project - you will work in groups of 2**. Sign up for a
team in the teammate signup sheet below. If you have a person you’d
prefer to work with, sign up as a team, otherwise just put your name on
the sheet – either as Team Member 1 on a new team, or as Member 2 next
to another solo student. If you can’t find a team by the end of the week
this is assigned, let our teaching team know  and we will help you find
one. We will be spending a significant in-class time doing "workshop"
sessions on the assignment with your teaching team present.

**Basic procedure**

1.  Choose City X. Choose a city with similar hydrology, topography, and
    climate, with available data comparable to Calgary's. Don't choose a
    city on the ocean, or a city in a desert.

2.  Gather open data from Calgary’s open data site
    (<https://data.calgary.ca/>) and City X's open data site as well as
    other internet sources. Build a fishnet grid for each city with the
    same sized cells. Calgary flood extents and other data are hosted on
    canvas in Module 5, they are called midTermProject\_Data.zip

3.  Using what we’ve learned about feature engineering, build as many
    useful variables describing the natural, hydrological and built
    environment factors that might help explain flood inundation. Form a
    hypothesis to motivate your search for data. **You must include at
    least one feature from the watershed analysis. **You need to build
    the exact same features for Calgary AND City X. We will have a
    feature engineering workshop in class. 

4.  Join your features to the vector Fishnet.

5.  Move your Fishnet datasets into R and set up your model workflow
    based on our Preservation markdown.

6.  Build logistic regressions for Calgary. Make **a training set and
    test your model by predicting on a test set**.  Set up a workflow
    with goodness of fit metrics - a confusion matrix for example.
    Experiment until you find a model with decent predictive
    performance.

7.  Using the predict function, feed *the exact same independent
    variables* from City X into your model object. Output a prediction
    for each cell in that city.

8.  Visualize and communicate your results.

**Deliverables**

A technical report - maximum of 4 pages - in pdf form (*or *knitted R
Markdown document html) due on Canvas.

**The audience for your report is the Director. Assume they are
technically savvy - they were the ones who asked you to do regression
modeling for this task!** Use clear language, good organization and
simple interpretation.

Your report should contain:

1.  An introduction framing the planning motivation for your algorithm
    and presenting the top-line results. How does this modeling approach
    work? Why are we doing it this way? What was the result of the
    analysis?

2.  A section outlining your feature selection. This should contain maps
    or plots related to *four *of your model features. Annotate with
    some commentary about why these are useful or un-useful
    features. **This must include at least one watershed feature (even
    if it’s not significant).**

3.  A section on your model. This should contain your final logistic
    regression model summary table, your ROC curve plot, and a table of
    confusion metrics (True Positives, False Positives, etc). Please
    discuss the modeling process, the training and testing process,
    comment on the accuracy of the model, and make a table that uses
    simple statements to describe the four model outcomes ("A True
    Positive is...").

4.  Support your modeling section with three maps:

    1.  One showing true positives, true negatives, false negatives and
        false positives for the **training set in Calgary.**

    2.  **S**econd, your inundation predictions for Calgary (entire
        dataset);

    3.  Third, predictions for your comparable city.

5.  A conclusion section commenting on the quality of your model, its
    usefulness in your comparable city, and your opinion about the
    utility of this modeling approach.

<img src="media/image1.png" style="width:5.76806in;height:2.98819in" />

There are three files as illustrated above in Module 5 on Canvas
(midtermdata.zip). There is **the satellite image with flood
inundation** in light blue and cloud cover in red. You will have to
reclassify into a raster of 0=no inundation and 1=inundation. **The pink
polygon is the study area** which you will note has a greater area than
the raster. You can 1) create a fishnet at the study area extent in
ArcGIS with [this toolLinks to an external
site.](https://desktop.arcgis.com/en/arcmap/10.3/tools/data-management-toolbox/create-fishnet.htm),
or in R with [sf::st\_make\_grid.Links to an external
site.](https://www.rdocumentation.org/packages/sf/versions/0.8-0/topics/st_make_grid) You
will then have to 2) create an unique id field in the fishnet shapefile
and 3) use zonal statistics as table to get the raster values into the
fishnet. Finally 4) remove the grid cells to the north for which we have
no raster data from the fishnet.  **You also have a DEM raster** as
above. For your comparison city, look to their open data site for
Digital Elevation Map (DEM). Otherwise, go to EarthExplorerLinks to an
external site. and download Shuttle Radar Topography Mission data
([SRTMLinks to an external
site.](https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-shuttle-radar-topography-mission-srtm-1-arc?qt-science_center_objects=0#qt-science_center_objects))
data.

**Project management**

You and your partner need to collect **the same features** (with the
same variable names) for both cities. To make a model work – your comp
city needs to have the *exact same features with the exact same variable
names so that the predict function can run using your model. *Keep in
mind that your model has specific parameters – if you can’t collect the
same data for City X as you can for Calgary (and for your model), it
won’t work!

I suggest you split your cities by teammate and create a shared table on
google docs for data tasks. This 'data matrix' will help you keep track
of what data is being collected by whom. 

A note about time management – if you have not done a project of this
nature before (in MUSA 508 for instance), you might not realize that it
is *very *hard to predict when and how things might go wrong with your
code. It is imperative that you start early.

**Assignment 3 Rubrics**

<table>
<tbody>
<tr class="odd">
<td><strong>Criteria</strong></td>
<td><strong>Ratings</strong></td>
<td><strong>Pts</strong></td>
</tr>
<tr class="even">
<td><p>This criterion is linked to a Learning OutcomeIntroduction and Planning Motivation</p>
<p>The Planning motivation for your algorithm and an overview of the modeling process and outcomes.</p></td>
<td></td>
<td>5 pts</td>
</tr>
<tr class="odd">
<td><p>This criterion is linked to a Learning OutcomeAnnotated Maps</p>
<p>Four of your more original (as you deem it), yet statistically significant features. Annotate as you see fit.</p></td>
<td></td>
<td>5 pts</td>
</tr>
<tr class="even">
<td><p>This criterion is linked to a Learning OutcomeModels &amp; Feature Engineering</p>
<p>Model outputs and basic goodness of fit with annotations, annotations and plots associated with feature engineering.</p></td>
<td></td>
<td>5 pts</td>
</tr>
<tr class="odd">
<td><p>This criterion is linked to a Learning OutcomeError Metrics &amp; Comp City</p>
<p>Maps of your predictions in training/testing and your "challenge" city</p></td>
<td></td>
<td>5 pts</td>
</tr>
<tr class="even">
<td><p>This criterion is linked to a Learning OutcomeCompleteness</p>
<p>Was it on time? Does it contain all of the key elements specified on the list of deliverables? Is it the specified length?</p></td>
<td></td>
<td>5 pts</td>
</tr>
</tbody>
</table>

Total Points: 25
