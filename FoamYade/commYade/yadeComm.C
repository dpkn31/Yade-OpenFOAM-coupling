#include "yadeComm.H" 

Foam::yadeComm::yadeComm()
{
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
}


void Foam::yadeComm::send_data(int recvRank, std::vector<double>& data) 
{
  MPI_Send(&data.front(), data.capacity(), MPI_DOUBLE, recvRank, sendTag, MPI_COMM_WORLD); 
}


void Foam::yadeComm::recv_data(int sendRank, std::vector<double>& data)
{

  MPI_Recv(&data.front(),data.capacity(),MPI_DOUBLE,sendRank,sendTag, MPI_COMM_WORLD,&status); 

}


void Foam::yadeComm::sendOneDouble(int recvRank, double& data)
{
  MPI_Send(&data,1,MPI_DOUBLE,recvRank,sendTag, MPI_COMM_WORLD); 
}



void Foam::yadeComm::cast_double_array_data(int cast_rank, std::vector<double>& array)
{
  MPI_Bcast(&array.front(),array.capacity(),MPI_DOUBLE, cast_rank, MPI_COMM_WORLD); 

}


void Foam::yadeComm::cast_one_double(int cast_rank,double& value) 
{ 
  MPI_Bcast(&value,1,MPI_DOUBLE, cast_rank, MPI_COMM_WORLD); 
}



void Foam::yadeComm::cast_integer_data(int cast_rank,int& value) 
{ 
  MPI_Bcast(&value,1,MPI_INT, cast_rank, MPI_COMM_WORLD); 
}

void Foam::yadeComm::procReduceMaxInt(int& sendValue, int& recValue)

{
  
  MPI_Allreduce(&sendValue,&recValue,1, MPI_INT, MPI_MAX,MPI_COMM_WORLD); 
  
}

void Foam::yadeComm::procReduceSumDouble(double& sendValue, double& recValue){
  MPI_Allreduce(&sendValue,&recValue,1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD); 
} 








