# Celltracker results analyzer.

Author: Guillermo Prol Castelo. Contact at gprolcast@gmail.com.

For MATLAB.

August 30th, 2019.

## Code description

When we are done tracking, we can save the progress in a *.mat* file, which has a spreadsheet-like structure with four untitled columns containing, from left to right, the X and Y positions, frame number and track ID (i.e. the first cell tracked has an ID value of 1, the second cell is 2, and so on).

This file can be opened on a MATLAB script and consequently analyzed. We begin by specifying the parameters: number of frames, time interval in seconds and the pixels to microns ratio. We also define the location path of the files we obtained above. Once that is done we create two matrices (one for X positions and another one for Y positions) such that cells positions are located on columns and the correspondent frame on a row. That is:

|         | cell 1 | cell 2 | cell 3 | etc. |
| :------ | ------ | ------ | ------ | ---- |
| frame 1 |        |        |        |      |
| frame 2 |        |        |        |      |
| frame 3 |        |        |        |      |
| etc     |        |        |        |      |

At the same time, a variable named *dummy\_frames* saves the number of frame that the first position of a cell corresponds to. This helps to then classify the cells into primary and daughter cells: the former are those whose first frame is less than or equal to 16, and the rest are daughter cells. The length in frames of the matrix containing all cells is obtained and a histogram of these is plotted. Furthermore, on the command window the average frame length is displayed (both in frames and seconds).

Velocities are then calculated for each of the positions matrices as the difference of positions on frame *i+1* and *i*, divided by the time interval specified at the beginning of the script. From these we can get the modulus of the velocity and, then, the average for each frame interval. These averages are plotted, one for each type of cells -- so there will be a total of three average velocity plots: for all cells, primary and daughter cells.

The MATLAB \*msdanalyzer* class is used to calculate and plot mean square displacements and velocity autocorrelations. Firstly, though, our data is structured so that *msdanalyzer* can understand it (i.e. in the MATLAB class *cell*). These calculations are done for all cells, primary and secondary cells, but the velocity autocorrelation is not working for the primary and secondary cells.

Finally, a chemotaxis plot is drawn with the output of a change of coordinates of our original matrices of positions. Each of a cell's positions is subtracted from their respective origin; that is if cell *i* begins at position
$$
x_{1i}
$$
then that position is now *0* and the rest are 
$$
x_{ji}-x_{1i}.
$$
We plot these positions with the help of  *msdanalyzer*.

## Issues

- The velocity autocorrelations calculation is not being done for the daughter cells. 
- There seems to be no distinction between the velocity autocorrelations for the primary cells and all of the cells.