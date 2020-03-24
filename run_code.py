import time
import os

time_start = time.time()


#remember to change the location of your boost directory. In my case, as in the example below, I have "/home/mat/Documents/Libraries/boost/"
os.system ("g++ -std=c++11 -O3 -o UrbanLikeEnvironment UrbanLikeEnvironment.cc  -L/home/mat/Documents/Libraries/boost/ -I/home/mat/Documents/Libraries/boost/")
os.system ("./UrbanLikeEnvironment ")


elapsed = time.time() - time_start
print ("cpu time: ", elapsed)
