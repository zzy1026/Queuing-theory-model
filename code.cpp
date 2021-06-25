#include <iostream>
#include <list>
#include <cmath>
#include <math.h>
#include <algorithm> 

using namespace std;


// randomly calculates negative-exponenetially-distributed-time

int random(long seed )
{   
    long a=seed;
    long b;
    b=(seed*(16807))%(2147483647);
    return b;
}

double nedt(double rate,long &seed)
{
    double u;
     
    seed=random(seed);
    u = (double)seed/2147483647;
    return ((-1/rate)*log(1-u));

}

// Mean and Standard Deviation 
double ave(double a[],int n)
{
    double sum=0;
    for (int i=0;i<n;i++)
    sum+=a[i];
    return sum/n;
}
double sd(double a[],int n)
{
    double sum=0,v,av;
    av=ave(a,n);
    for (int i=0;i<n;i++)
    sum=(a[i]-av)*(a[i]-av)+sum;
    v=sum/n;
    return sqrt(v);
} 



class Event {
	double eventTime;
	bool isArrival; 
    int num;

public:
	Event(double etime, bool arrival, int number) {
		eventTime = etime;
        isArrival = arrival;
        num=number;
	}

	double getEventTime() {
		return eventTime;
	}
        
    int getNum() {
        return num;
    }
	bool getIsArrival() {
		return isArrival;
	}

	bool operator==(const Event &rhs) const {
        return rhs.eventTime == eventTime;
    }

	bool operator>=(const Event &rhs) const {
        return rhs.eventTime >= eventTime;
    }

    bool operator>(const Event &rhs) const {
        return rhs.eventTime > eventTime;
    }
};


class GEL { // Global Event List

	std::list<Event> GlobalEventList;

public:
	GEL() {
        GlobalEventList = std::list<Event>();
    }
	void insert(Event event) {
        //cout << "begin insert" << endl;

        //cout << "insert event.isArrival: " << event.getIsArrival() << endl;
        //cout << "size = " << GlobalEventList.size() << endl;

		if (GlobalEventList.size() == 0) {
            GlobalEventList.push_front(event);
            return;
        }

		for (std::list<Event>::iterator itr = GlobalEventList.begin(); itr != GlobalEventList.end(); itr++) {
			if (itr->getEventTime() > event.getEventTime()) {
                GlobalEventList.insert(itr, event);
                return;
            }
        }

        GlobalEventList.push_back(event);
	

        //cout << "end insert " << endl;

	} // insert sorted by events time

	Event removeFirst() {

        //cout << "begin removeFirst" << endl;
		Event firstElement = GlobalEventList.front();
		GlobalEventList.pop_front();

        //cout << "end removeFirst" << endl;

		return firstElement;
	}
};

int  main()
{
    // should be read in from command line
    double lambda;
    double mu;
    double maxbuffer;
    long seed1=1814472367,seed2=1328130192;

    cout << "lambda: ";
    cin >> lambda;

    cout << "mu: ";
    cin >> mu;

    cout << "Buffer Size: ";
    cin >> maxbuffer;

    // initialize
    int length = 0; //queue length
    int dropNum = 0;
    double sumLength = 0;
    double time = 0;
    double busy = 0;//work time of servers
    double serviceTime = 0; //service time
    int N=30000; //event number
    double arr[N]; //arrival time
    double dep[N]; //departure time
    int server=2; //server number
    int count=1000;//customer number
    int k=0,m=0;

    //first arrival
    Event e = Event(time + nedt(lambda,seed1), true,1);
    GEL eventList = GEL();
    eventList.insert(e);

    // process event
    // just going by the number given

    for (int i = 0; i < N; i++)
    {
        //if(server!=0 && server!=1 && server!=2)
        //cout<<"error   "<<server<<endl;
        // get closest event and update time
        e = eventList.removeFirst();

        // sums length by multiplying length by elapsed time
        sumLength += max(0, length - 1) * (e.getEventTime() - time);
        //cout << "prev time: " << time  << "   event time: " << e.getEventTime() << endl; 

        // updates time
        time = e.getEventTime();

        // handles Arrival event
        if (e.getIsArrival())
        {
            

            //cout << "is Arrival, i: " << i << endl;
            // insert new arrival event
            eventList.insert(Event(time + nedt(lambda,seed1), true, k+1));

            //cout << "length: " << length << endl;

            // if server is free, schedule a departure event, and update length
            if (length == 0)
            {
            arr[k++]=e.getEventTime();
            //cout<<k<<"arr"<<'\t'<<e.getEventTime()<<'\n';
            if (server == 2) //there are two servers free
            {
                //cout << "hello from length = 0" << endl;
                serviceTime = nedt(mu,seed2);
                //cout << "service time: " << serviceTime << endl;
                busy += serviceTime;
                dep[m++]=time + serviceTime;
                //cout<<m<<"dep"<<'\t'<<dep[m-1]<<'\n';
                eventList.insert(Event(time + serviceTime, false ,m));
                server--;

            }

            // if server is free, schedule a departure event, and update length
            else if (server == 1)//there are only one servers free
            {
                //cout << "hello from length = 0" << endl;
                serviceTime = nedt(mu,seed2);
                //cout << "service time: " << serviceTime << endl;
                busy += serviceTime;
                dep[m++]=time + serviceTime;
                //cout<<m<<"dep"<<'\t'<<dep[m-1]<<'\n';
                eventList.insert(Event(time + serviceTime, false, m));
                length ++;
                server--;
                // this assumes maxbuffer is at least one, 
                // which is a good assumption because no buffer 
                // would have max buffer equal to 1
            }
            }
            // else if room in buffer
            // maxbuffer = -1 denotes infinite buffer
            else if (maxbuffer == -1 || length - 1 < maxbuffer)
            {
            arr[k++]=e.getEventTime();
            //cout<<k<<"arr"<<'\t'<<e.getEventTime()<<'\n';
                length ++;
                // handles generating of service time time when departure event created
            }
            else // no room in buffer
            {
                dropNum ++;
            }
        }

        // handles departure event
        else 
        {   
            //cout<<e.getNum()<<"dep"<<'\t'<<dep[e.getNum()-1]<<'\n';
            //cout << "is departure" << endl;
            server ++;
            if(server!=2) length --;
            // if customer still in queue, create a departure event
            if (length > 0)
            {
                serviceTime = nedt(mu,seed2);
                //cout << "service time: " << serviceTime << endl;

                busy += serviceTime;
                dep[m++]=time + serviceTime;
                //cout<<m<<"dep"<<'\t'<<dep[m-1]<<'\n';
                eventList.insert(Event(time + serviceTime, false ,m));
                server--;
            }
        }
    if(m==count) break; 

    }
    
    count= min(m,k);
    cout<<"customer number:"<<count<<'\n';

    // calculate system time
    double sys[count]; //system time
    for(int j=0;j<count;j++)
    { 
        sys[j]=dep[j]-arr[j];
    }
    
    double mean=ave(sys,count);
    cout<<"Mean system time:"<<mean<<endl;
    double s=sd(sys,count);
    cout<<"Standard Deviation of system time:"<<s<<endl;
    double h=1.962*s/sqrt(count);
    cout<<"the 95% confidence interval:("<<mean-h<<","<<mean+h<<")"<<endl;
    
    if(busy<(2*time)) cout<<"Utilization: " << busy / (2*time) << endl;
    else cout << "Utilization: " << 1.0f << endl;
    cout << "Mean queue length: " << sumLength / time << endl;
    cout << "Number of  dropped: " << dropNum << endl << endl << endl;


	return 0;
}

