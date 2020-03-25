def defineheader(TBU):
# defines information for the header of a ToF-GUI readable file
  count=TBU.shape[0] #count of (T,B,U) values
  headerstring = ''      # write header first to string, then at the end of outer loop to file

  outputformat = '%s\t';
  for i in range(count):
    outputformat = outputformat + ' %3gK\t'
  outputformat = outputformat[0:-1] + '\n'
  headerstring = headerstring + outputformat % (('SetTemp',)+tuple(TBU[:,0]))
  
  outputformat = '%s\t';
  for i in range(count):
    outputformat = outputformat + ' %5.5gK\t'
  outputformat = outputformat[0:-1] + '\n'
  headerstring = headerstring + outputformat % (('RealTemp',)+tuple(TBU[:,0]))
  
  outputformat = '%s\t';
  for i in range(count):
    outputformat = outputformat + ' %5.3gT\t'
  outputformat = outputformat[0:-1] + '\n'
  headerstring = headerstring + outputformat % (('SetMag',)+tuple(TBU[:,1]))
  
  outputformat = '%s\t';
  for i in range(count):
    outputformat = outputformat + ' %5.3gT\t'
  outputformat = outputformat[0:-1] + '\n'
  headerstring = headerstring + outputformat % (('RealMag',)+tuple(TBU[:,1]))
  
  outputformat = '%s\t';
  for i in range(count):
    outputformat = outputformat + ' %3gV\t'
  outputformat = outputformat[0:-1] + '\n'
  headerstring = headerstring + outputformat % (('SetBias',)+tuple(TBU[:,2]))
  
  outputformat = '%s\t';
  for i in range(count):
    outputformat = outputformat + ' %6.6gV\t'
  outputformat = outputformat[0:-1] + '\n'
  headerstring = headerstring + outputformat % (('SetMag',)+tuple(TBU[:,2]))
  
  # outputformat for datamatrix
  outputformat = '%5.3e\t';
  for i in range(count):
    outputformat = outputformat + ' %12.10f\t'
  outputformat = outputformat[0:-1] + '\n'
 
  return outputformat, headerstring

