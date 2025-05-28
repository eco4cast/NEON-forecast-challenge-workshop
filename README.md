[![DOI](https://jose.theoj.org/papers/10.21105/jose.00259/status.svg)](https://doi.org/10.21105/jose.00259) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14018090.svg)](https://doi.org/10.5281/zenodo.14018090)


# NEON Forecast Challenge Workshop Pre-Workshop Instructions

This is a repository for materials to complete the NEON Ecological Forecasting Challenge tutorial. This tutorial repository is a dynamic repository and may change in future iterations.

The tutorial is designed as a 90-minute, standalone session that includes an introductory presentation, a guided demonstration of forecast code, and open time for discussion. If a longer time period is available, this repository has additional activities and content that includes more advanced topics. These materials can also be completed as a self-paced tutorial.

## 1. Pre-workshop intructions

### Setting up your R environment

R version 4.2 is required to run the code in this workshop. You should also check that your Rtools is up to date and compatible with R 4.2, see (<https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html>). Although instructions assume users have Rstudio installed, this is not a strict requirement.

The following packages need to be installed using the following code.

```{r}
install.packages('remotes')
install.packages('tidyverse') # collection of R packages for data manipulation, analysis, and visualisation
install.packages('lubridate') # working with dates and times
remotes::install_github('eco4cast/neon4cast') # package from NEON4cast challenge organisers to assist with forecast building and submission
```

### Get the code

There are 3 options for getting the code locally so that you can run it, depending on your experience with Github/Git you can do one of the following:

1.  **Fork (recommended)** the repository to your Github and then clone the repository from your Github repository to a local RStudio project. This will allow you to modify the scripts and push it back to your Github.

-   Find the fork button in the top right of the webpage --\> Create Fork. This will generate a copy of this repository in your Github.
-   Then use the \<\> Code button to copy the HTTPS link (from your Github!).
-   In RStudio, go to New Project --\> Version Control --\> Git.
-   Paste the HTTPS link in the Repository URL space, and choose a suitable location for your local repository --\> Create Project.
-   Open the .Rmd file

2.  **Clone** the workshop repository to a local RStudio project. Your local workspace will be set up and you can commit changes locally but they won't be pushed back to the Github repository.

-   Use the \<\> Code button to copy the HTTPS link.
-   In RStudio go to New Project --\> Version Control --\> Git.
-   Paste the HTTPS link in the Repository URL space, and choose a suitable location for your local repository --\> Create Project.
-   Open the .Rmd file

3.  **Download** the zip file of the repository code. You can save changes (without version control) locally.

-   Find the \<\> Code button --\> Download ZIP.
-   Unzip this to a location on your PC and open the `NEON-forecast-challenge-workshop.Rproj` file in RStudio.

More information on forking and cloning in R can be found at [happygitwithr](https://happygitwithr.com/fork-and-clone.html), a great resource to get you started using version control with RStudio.

For the workshop you can follow along via the rmarkdown document (`Submit_forecast/NEON_forecast_challenge_workshop_aquatics.Rmd`) or the md (`Submit_forecast/NEON_forecast_challenge_workshop_aquatics.md`), both of which can be downloaded here or you can fork the whole repository.

## 2. Introductory presentation

An introductory presentation that covers key forecasting concepts and a broad introduction the the Challenge and NEON data is provided with this repository as both a [PowerPoint](Introductory_presentation.pptx) and a [PDF](Introductory_presentation.pptx). This presentation can be tailored to the audience based on expected familiarity with forecasting concepts and the data being forecasted.

## 3. Coding walk-through

The hands-on materials are split into three sections, each in a separate subfolder:

1.  [`Submit_forecast`](Submit_forecast) - want to get started with making and submitting forecasts to the Challenge - start here!
2.  [`Automate_forecasts`](Automate_forecasts) - you've made a model and submitted a forecast successfully? Automate your workflow to submit a new forecast every day
3.  [`Analyse_scores`](Automate_forecasts) - interested in knowing how different forecasts are performing within the Challenge? This is a brief introduction in how to access the scored forecasts and do some simple analyses.

For a 90-minute workshop, the hands-on R markdown walk through uses the materials in the `Submit_forecast` directory (30-40 minutes). The Rmarkdown includes descriptive text to reinforce some of the ideas introduced in the presentation and guides participants through the stages of forecast development from importing 'target' data (observations), plotting and investigating data, developing a model, generating a forecast, adjusting the format to Challenge standards, and submitting their forecast to the Challenge.

Materials in the other directories are recommended for later completion and cover more advanced topics - see [Tutorial extension](#5-tutorial-extension) section below.

## 4. Open time

For a standard 90-minute workshop, time at the end of the session can be used to debug user-specific code issues (package installation etc.), allow participants to modify the forecast model code available in the template ([Submit_forecast/forecast_code_template.R](Submit_forecast/forecast_code_template.R), explore alternative forecast methods, or form forecasting teams to submit additional forecasts to the Challenge beyond the workshop session. Some potential topics and additional activities are detailed at the end of the Rmarkdown, including code for alternative forecasting methods or possible modifications to the example forecast.

## 5. Tutorial extension

To extend the tutorial beyond 90 minutes, we have provided materials that walk participants through how forecast submission could be automated (`Automate_Forecast`). Automation is a crucial step for iterative, near-term ecological forecasting, and is highly encouraged when participating in the Challenge after completing the tutorial. Finally, materials are available that introduce participants to how a forecast is evaluated and how the scores (evaluated forecast performance) can be analysed (`Analyse_Scores`). These materials could also be used in a self-paced manner after a workshop.

------------------------------------------------------------------------

### Optional: Use Docker

An alternative to running the materials locally, is to use a Docker container that has all the packages pre-installed.

#### Installing Docker

Go to <https://docs.docker.com/get-docker/> to install the relevant install for your platform (available for PC, Mac and Linux). Also see <https://docs.docker.com/desktop/>.

NOTE: \* If you're running Windows, you will need WSL (Windows Subsystem for Linux) \* If you're running a Linux distribution, you may have to enable Viritualization on your computer (see [here](https://stackoverflow.com/questions/76646465/unable-to-launch-docker-desktop-on-ubuntu/76655270#76655270))

#### Running a docker container

1.  Launch Docker Desktop (either from the Command Line or by starting the GUI)
2.  At the command line run the following command which tells docker to `run` the container with the name `eco4cast/rocker-neon4cast` that has all the packages and libraries installed already. The `PASSWORD=yourpassword` sets a simple password that you will use to open the container. The `-ti` option starts both a terminal and an interactive session.

```         
docker run --rm -ti -e PASSWORD=yourpassword -p 8787:8787 eco4cast/rocker-neon4cast
```

This can take a few minutes to download and install. It will be quicker the next time you launch it.

3.  Open up a web browser and navigate to `http://localhost:8787/`
4.  Enter the username: `rstudio` and password: `yourpassword`
5.  You should see a R Studio interface with all the packages etc. pre-installed and ready to go.

You can close this localhost window (and then come back to it) but if you close the container from Docker (turn off your computer etc.) any changes will be lost unless you push them to Github or exported to your local environment.

Then follow the instructions for Getting the Code above (`2. Get the Code`).

#### Note for Mac users

If you are running Docker on ARM architecture (e.g., Apple Silicon M chips), the following Docker command is needed

```         
docker run --platform linux/amd64 --rm -ti -e PASSWORD=yourpassword -p 8787:8787 eco4cast/rocker-neon4cast
```

If you are using macOS Sequoia (and future macOS versions), you may get this error:

```         
rosetta error: Rosetta is only intended to run on Apple Silicon with a macOS host using Virtualization.framework with Rosetta mode enabled
```

If you get this error, you will need to install `rosetta` before running the docker command above.

```         
softwareupdate --install-rosetta
```
