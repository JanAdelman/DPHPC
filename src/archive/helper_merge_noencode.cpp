#include <iostream>
#include <cstring>
#include <algorithm>
#include <string>
#include <vector>
#include <tuple>
#include <type_traits>
#include <ostream>
#include <fstream>

#define K_SIZE 2

struct tuple_t
{
    int idx;
    char seq[K_SIZE];
};


struct triple_t
{
    int idx;
    int b;
    int b2;
};

struct tuple_ISA{
    int SA;
    int B;
};

int bucket_id(int* displ, tuple_ISA &SA_B, int world_size){
    for(int i=0;i<world_size-1;i++){
        if (displ[i]<=SA_B.SA && displ[i+1]>SA_B.SA){
            //std::cout<<"SA "<<SA_B.SA<<"i "<<i<<std::endl;
            return i;
        }
    }
    return world_size-1;
}

int bucket_id_shift(int* displ, int B, int world_size, int* sendcounts){
    if(B>displ[world_size-1]+sendcounts[world_size-1]-1){
        return -1;
    }
    for(int i=0;i<world_size-1;i++){
        if (displ[i]<=B && displ[i+1]>B){
            //std::cout<<"SA "<<SA_B.SA<<"i "<<i<<std::endl;
            return i;
        }
    }
    return world_size-1;
}

void rebucketing(std::vector<tuple_ISA> &SA_B, std::vector<tuple_t> &kmers, size_t size, int* displ, int* counts, int world_rank, int world_size){
    SA_B[0].B=displ[world_rank];
    SA_B[0].SA=kmers[0].idx;
    int id_zero=bucket_id(displ,SA_B[0],world_size);
    counts[id_zero]+=1;
    for(int i=1;i<size;i++){
        if(!std::lexicographical_compare(kmers[i-1].seq, kmers[i-1].seq + K_SIZE, kmers[i].seq, kmers[i].seq + K_SIZE)){
            SA_B[i].B=SA_B[i-1].B;
            SA_B[i].SA=kmers[i].idx;
            int id=bucket_id(displ,SA_B[i],world_size);
            counts[id]+=1;
        }
        else{
            SA_B[i].B=displ[world_rank]+i;
            SA_B[i].SA=kmers[i].idx;
            int id=bucket_id(displ,SA_B[i],world_size);
            counts[id]+=1;
        }
    }
}

//TRIPLET
void rebucketing(std::vector<tuple_ISA> &SA_B,std::vector<triple_t> &input, size_t size, int* displ, int world_rank, int world_size,int* counts){
    SA_B[0].B=displ[world_rank];
    SA_B[0].SA=input[0].idx;
    int id_zero=bucket_id(displ,SA_B[0],world_size);
    counts[id_zero]+=1;
    for(int i=1;i<size;i++){
        if((input[i-1].b == input[i].b and input[i-1].b2 == input[i].b2)){
            SA_B[i].B=SA_B[i-1].B;
            SA_B[i].SA=input[i].idx;
            int id=bucket_id(displ,SA_B[i],world_size);
            counts[id]+=1;
        }
        else{
            SA_B[i].B=displ[world_rank]+i;
            SA_B[i].SA=input[i].idx;
            int id=bucket_id(displ,SA_B[i],world_size);
            counts[id]+=1;
        }
    }
}


std::vector<tuple_ISA> probing_alltoallv(std::vector<tuple_ISA> &sendbuf, int* sdispls, int size, int world_size, int* sendcounts, MPI_Comm comm, int world_rank,MPI_Datatype TYPE){


    //Send to every neighbour
    for (int i = 0; i < world_size; i++){
        if(i!=world_rank){
            //std::cout << "Sending out" << world_rank <<" " << sendcounts[i]<<std::endl;
            MPI_Request request;
            MPI_Isend(&sendbuf[sdispls[i]],sendcounts[i], TYPE, i, 0, comm, &request);
            //std::cout << world_rank << " sent in probing"<<std::endl;
        }
    }


    //Allocate array of global size on this processor
    std::vector<tuple_ISA> recv_buffer(size);
    //Size of recieve buffer used so far
    int used_size = 0;
    int recieved_size;
    //Recieve from every neighbour
    for (int i = 0; i < world_size; i++){
        if(i!=world_rank){
            MPI_Status status;
            MPI_Recv(&recv_buffer[used_size], size, TYPE, i, 0, comm, &status);
            //MPI_Sendrecv(&sendbuf[sdispls[i]], sendcounts[i], TYPE, i, 0, );
            MPI_Get_count(&status, TYPE, &recieved_size);
            //Increase used by recieved size to prepare next offset
        }
        else{
            std::copy(&sendbuf[sdispls[world_rank]],&sendbuf[sdispls[world_rank]+sendcounts[world_rank]],&recv_buffer[used_size]);
            recieved_size=sendcounts[world_rank];
        }
        //Increase used by recieved size to prepare next offset
        used_size += recieved_size;

    }

    return recv_buffer;
}



void print_char_array(char *input, size_t size)
{
    for (int i = 0; i < size; ++i)
        std::cout <<"   "<< input[i];
    std::cout << std::endl;
}

void print_int_array(std::vector<int> &input, size_t size)
{
    for (int i = 0; i < size; ++i)
        std::cout << input[i]<<",";
    std::cout << std::endl;
}

void get_kmers_adapt(const char *input, const int k, std::vector<tuple_t> &kmers, size_t size, int displacement)
{
    for (int i = 0; i < size; i++)
    {
        memcpy(kmers[i].seq, input + i, k);
        kmers[i].idx = displacement+ i;
    }
}

bool char_array_comp(const char* a, const char* b, int size){
    bool same = true;
    for (int i = 0; i < size; i++) {
        if (*(a + i) != *(b + i))
            same = false;
    }
    return same;
}

bool tuple_t_compare(const tuple_t &a, const tuple_t &b) {
    int res = std::strcmp(a.seq, b.seq);
    if (res < 0)
        return true;
    else if (res == 0)
        return a.idx<b.idx;
    else 
        return false;
}


bool triple_t_compare(const triple_t &a, const triple_t &b)
{
    /*
    if (a.b < b.b)
        return true;
    else if (a.b > b.b)
        return false;
    else
    {
        if (a.b2 < b.b2)
            return true;
        else if (a.b2 > b.b2)
            return false;
        else
        {
            if (a.idx < b.idx)
                return true;
            else // if(std::get<2>(a)>std::get<2>(b))
                return false;
        }
    }
    */
   return a.b < b.b || (a.b == b.b && a.b2 < b.b2) || (a.b == b.b && a.b2 == b.b2 && a.idx < b.idx);
}



void t_sort(std::vector<tuple_t> &input, size_t size)
{
    std::sort(input.begin(), input.end(), tuple_t_compare);
}
void t_sort(std::vector<triple_t> &input, size_t size)
{
    std::sort(input.begin(), input.end(), triple_t_compare);
}

void t_print(std::vector<tuple_t> &input, size_t size)
{
    for (int i = 0; i < size; i++){
        print_char_array(input[i].seq, K_SIZE);
        std::cout<<"     "<<input[i].idx<<std::endl;
    }
    std::cout << "---" << std::endl;
}
void t_print(std::vector<triple_t> &input, size_t size)
{
    for (int i = 0; i < size; i++)
    {
        std::cout << "    B: " << input[i].b << std::endl;
        std::cout << "    B2: " << input[i].b2 << std::endl;
        std::cout << "    SA: " << input[i].idx << std::endl;
    }
    std::cout << "---" << std::endl;
}

void t_print_flat(const triple_t *input, size_t size, int rank, int target)
{
    std::cout <<rank<<"(from):(to)"<<target<< "->";
    for (int i = 0; i < size; i++)
    {
        std::cout << input[i].b <<","<< input[i].b2<<","<<input[i].idx<<";";//  << std::endl;
    }
    std::cout << std::endl;
}


void t_print(std::vector<tuple_ISA> input, size_t size)
{
    for (int i = 0; i < size; i++)
    {
        std::cout << " B: " << (input.begin() + i)->B << std::endl;
        std::cout << " SA: " << (input.begin() + i)->SA << std::endl;
    }
    std::cout << "---" << std::endl;
}

void debug_tuple_print(std::vector<tuple_ISA> &input, size_t size)
{

    std::ofstream outfile ("./result.txt");

    for (int i = 0; i < size; i++)
    {
        //std::cout << ((input + i)->SA) << ",";
        outfile << input[i].SA << ",";
    }

    outfile.close();
}

void reorder_to_stringorder(std::vector<int> &B,std::vector<tuple_ISA> &input,size_t size_inp, int displacement){
    for(int i=0;i<size_inp;i++){
        B[(input[i].SA)-displacement] = input[i].B;
    }
}



// Adapted from http://selkie-macalester.org/csinparallel/modules/MPIProgramming/build/html/mergeSort/mergeSort.html
std::vector<tuple_t> typename_t_sort(int height, int id, std::vector<tuple_t> &localArray, int size, MPI_Comm comm)
{

    MPI_Datatype MPI_TUPLE_STRUCT;
    int lengths[2] = {1, K_SIZE};
    MPI_Aint displacements[2] = {0, sizeof(int)};
    MPI_Datatype types[2] = {MPI_INT, MPI_CHAR};
    MPI_Type_create_struct(2, lengths, displacements, types, &MPI_TUPLE_STRUCT);
    MPI_Type_commit(&MPI_TUPLE_STRUCT);

    int parent, rightChild, local_height;
    std::vector<tuple_t> half1, half2, mergeResult;

    local_height = 0;

    half1 = localArray; // assign half1 to localArray

    while (local_height < height)
    { // not yet at top
        //std::cout << "Getting Chil " << id << std::endl;
        parent = (id & (~(1 << local_height)));
        //std::cout << "Got child " << id << std::endl;

        if (parent == id)
        { // left child

            rightChild = (id | (1 << local_height));

            //std::cout<<"height now"<<local_height<<"current id"<<id<<"right child"<<rightChild<<std::endl;

            int recieved_size;
            MPI_Status status;

            MPI_Probe(rightChild, 0, MPI_COMM_WORLD, &status);

            MPI_Get_count(&status, MPI_TUPLE_STRUCT, &recieved_size);
            // allocate memory and receive array of right child
            half2.resize(recieved_size);
            //std::cout<< "child from which "<<id<< "receives"<<rightChild<<std::endl;
            MPI_Recv(&half2[0], recieved_size, MPI_TUPLE_STRUCT, rightChild, 0, MPI_COMM_WORLD, &status);

           // mergeResult = (tuple_t *)malloc((size + recieved_size) * sizeof(tuple_t));
            // merge half1 and half2 into mergeResult

            mergeResult.resize(size+recieved_size);
            std::merge(half1.begin(), half1.end(), half2.begin(), half2.end(), mergeResult.begin(), tuple_t_compare);
            // reassign half1 to merge result
            half1 = mergeResult;

            if (local_height == height-1 && id == 0)
            {
                return half1;
            }

            size = size + recieved_size;
            //std::cout<<"current height: "<<local_height<<"max height"<<height<<std::endl;
            local_height++;

        }
        else
        { // right child
            // send local array to parent

            //std::cout << "Sedning " << id << std::endl;
            //std::cout<< "parent to which "<<id<< "send:"<<parent<<std::endl;
            MPI_Send(&half1[0], size, MPI_TUPLE_STRUCT, parent, 0, MPI_COMM_WORLD);
            local_height = height;
            //std::cout << "Sent " << id << std::endl;
        }
    }
    std::vector<tuple_t> empty (0);
    return empty;
}

void create_triple(std::vector<int> &B, std::vector<int> &B2, int size,int displ, std::vector<triple_t> &triple_arr){
    //triple_t triple_arr[size];
    for(int i=0;i<size;i++){
        triple_arr[i].b=B[i];
        //std::cout<<"B:"<<triple_arr[i].b<<std::endl;
        triple_arr[i].b2=B2[i];
        //std::cout<<"B2:"<<triple_arr[i].b2<<std::endl;
        triple_arr[i].idx=i+displ;
        //std::cout<<"SA:"<<triple_arr[i].idx<<std::endl;
    }
    //return triple_arr;

}
std::vector<triple_t> typename_t_sort(int height, int id, std::vector<triple_t> &localArray, int size, MPI_Comm comm)
{
    MPI_Datatype MPI_TRIPLE_STRUCT;
    int lengths_triple[3] = {1, 1, 1};
    MPI_Aint displacements_triple[3] = {0, sizeof(int), 2*sizeof(int)};
    MPI_Datatype types_triple[3] = {MPI_INT, MPI_INT, MPI_INT};
    MPI_Type_create_struct(3, lengths_triple,
                           displacements_triple, types_triple, &MPI_TRIPLE_STRUCT);
    MPI_Type_commit(&MPI_TRIPLE_STRUCT);

    int parent, rightChild, local_height;
    std::vector<triple_t> half1, half2, mergeResult;

    local_height = 0;

    half1 = localArray; // assign half1 to localArray

    while (local_height < height)
    { // not yet at top
        parent = (id & (~(1 << local_height)));

        if (parent == id)
        { // left child
            rightChild = (id | (1 << local_height));

            int recieved_size;
            MPI_Status status;
            MPI_Probe(rightChild, 0, comm, &status);

            MPI_Get_count(&status, MPI_TRIPLE_STRUCT, &recieved_size);
            // allocate memory and receive array of right child
            //half2 = (triple_t *)malloc(recieved_size * sizeof(triple_t));
            half2.resize(recieved_size);

            MPI_Recv(&half2[0], recieved_size, MPI_TRIPLE_STRUCT, rightChild, 0, comm, &status);
            //std::cout<<id<<rightChild<<std::endl;

            //t_print_flat(half2, recieved_size , rightChild, id);

            //mergeResult = (triple_t *)malloc((size + recieved_size) * sizeof(triple_t));
            mergeResult.resize(size+recieved_size);
            // merge half1 and half2 into mergeResult

            std::merge(half1.begin(), half1.begin() + size, half2.begin(), half2.begin() + recieved_size, mergeResult.begin(), triple_t_compare);

            // reassign half1 to merge result
            half1 = mergeResult;

            if (local_height == height-1 && id == 0)
            {
                return half1;
            }

            size = size + recieved_size;
            local_height++;
        }
        else
        { // right child
            // send local array to parent

            MPI_Send(&half1[0], size, MPI_TRIPLE_STRUCT, parent, 0, comm);

            local_height = height;
        }
    }
    std::vector<triple_t> empty (0);
    return empty;
}


    void naive_shift(std::vector<int> &input, const int h, MPI_Comm comm, const int world_rank, const int world_size, int *offsets_start, int local_length, int *offsets_end){
             //print_int_array(input, local_length);

            if (world_rank == 0)
            {
                //std::cout <<h<<" "<<local_length << std::endl;
                if (h < local_length){
                    std::copy(&input[h], &input[local_length], &input[0]);
                    //std::cout <<"Copying on " << world_rank << " "<< input[0] <<std::endl;
                }
                //print_int_array(input, local_length);
            }
            else{
                int shift_start = offsets_start[world_rank] - h;
                int shift_end = offsets_end[world_rank] - h;

                for (int rank = world_rank; rank >= 0; rank--){//Look at previous buckets

                    //std::cout<<world_rank <<" "<< shift_start<<" "<< shift_end << " "<< offsets_start[rank] <<" " <<offsets_end[rank] << std::endl;
                    if ((shift_start >= offsets_start[rank]) and (shift_start <= offsets_end[rank])){//Start is within the offsets of this rank
                        /*
                        std::cout << shift_end <<
                        if (shift_end >= offsets_start[world_rank]){//sending to the direct neighbour
                            MPI_Send(input, offsets_end[rank] - shift_start + 1, MPI_INT, rank, 0, comm);
                            std::cout << world_rank << "SENDING TO (1)" << rank << std::endl;
                            if (h < local_length)
                                std::copy(input + h, input + local_length, input);

                            break;
                        }
                        */
                        if (shift_end >= offsets_start[world_rank]){ //End is also within bounds -> Only send to one porcessor

                            //std::cout << world_rank << "SENDING TO (MODE 1)" << rank << std::endl;
                            //std::cout << "Sending size: " << local_length << " Inserting at: " << shift_start << std::endl;

                            MPI_Send(&input[0], h, MPI_INT, rank, shift_start, comm);

                            std::copy(&input[h], &input[local_length], &input[0]);
                            break;
                        }
                        else if (shift_end <= offsets_end[rank])
                        {
                            //std::cout << world_rank << "SENDING TO (MODE 2)" << rank << std::endl;
                            //std::cout << "Sending size: " << local_length << " Inserting at: " << shift_start << std::endl;

                            MPI_Send(&input[0], local_length, MPI_INT, rank, shift_start, comm);
                        }
                        else if (shift_end >= offsets_start[rank + 1])
                        {
                            //std::cout << world_rank << "SENDING TO (MODE 3)" << rank << std::endl;
                            //std::cout << "Sending size: " << local_length << " Inserting at: " << shift_start << std::endl;

                            MPI_Send(&input[0], offsets_end[rank] - shift_start + 1, MPI_INT, rank, shift_start, comm); //From calculated start to end of bucket -> Overlapping to next
                            MPI_Send(&input[offsets_end[rank] - shift_start + 1]
                                ,local_length - (offsets_end[rank] - shift_start + 1), MPI_INT, rank + 1, offsets_start[rank+1], comm); //MPI send end piece to rank where end is found; from start of bucket to calculated end
                            break;

                        }
                        else
                        {
                            std::cout << "Error 101 on rank "<<world_rank << std::endl;
                        }


                    }
                }
                if ((shift_start < 0) and (shift_end)>= 0)
                {
                        //std::cout << world_rank << "SENDING TO (MODE 4)" << 0 << std::endl;
                        MPI_Send(&input[local_length - (shift_end + 1)], shift_end + 1, MPI_INT, 0, 0, comm);
                }


            }


            //Recieving
            if (((offsets_end[world_size - 1] - offsets_end[world_rank]) >= h)) {//No zeroes are needed

                if (local_length > h){//Local shifting from direct neghbour
                    //MPI_Recv(input + local_length - h, h, MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//Recieve start MPI_recv(input,number_amount, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status)//Recieve start

                    //std::cout << "Awaiting recieve " <<world_rank<< std::endl;
                    MPI_Status status;
                    int position;
                    int number_amount;
                    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    position = status.MPI_TAG;
                    MPI_Get_count(&status, MPI_INT, &number_amount);


                    MPI_Recv(&input[position-offsets_start[world_rank]],number_amount, MPI_INT, MPI_ANY_SOURCE, position, MPI_COMM_WORLD, &status);//Recieve start
                    //std::cout << "Recieved " <<world_rank<< std::endl;
                }
                else
                {
                    //std::cout << "Awaiting recieve " <<world_rank<< std::endl;
                    MPI_Status status;
                    int position;
                    int number_amount;
                    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    position = status.MPI_TAG;
                    MPI_Get_count(&status, MPI_INT, &number_amount);

                    //std::cout << "position1" << position << std::endl;

                    MPI_Recv(&input[position-offsets_start[world_rank]],number_amount, MPI_INT, MPI_ANY_SOURCE, position, MPI_COMM_WORLD, &status);//Recieve start
                    //std::cout << "Recieved " <<world_rank<< std::endl;

                    if (number_amount < local_length){//If not enough recieve more

                        //std::cout << "Awaiting recieve more " <<world_rank<< std::endl;
                        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                        position = status.MPI_TAG;
                        MPI_Get_count(&status, MPI_INT, &number_amount);

                        //std::cout << "position2" << position << std::endl;

                        MPI_Recv(&input[position-offsets_start[world_rank]],number_amount, MPI_INT, MPI_ANY_SOURCE, position, MPI_COMM_WORLD, &status);//Recieve start
                        //std::cout << "Recieved more" <<world_rank<< std::endl;
                    }
                    /*
                    std::cout << "Recieving locally"<<world_rank<< std::endl;
                    MPI_Status status;
                    int number_amount;
                    MPI_Recv(input,local_length, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);//Recieve start
                    MPI_Get_count(&status, MPI_INT, &number_amount);
                    std::cout << "Amount " << world_rank <<" "<< number_amount <<std::endl;
                    if (number_amount < local_length){//If not enough recieve more

                        std::cout << world_rank << " waiting for more!" << std::endl;
                        MPI_Recv(input + number_amount,local_length - number_amount, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);//On a difffernt tagstream
                    }
                    */

                }
            }
            else if ((offsets_end[world_size -1] + 1 - offsets_start[world_rank] <= h))//Only zeroes are filled
            {
                //std::cout << "Zeroes are added to rank " << world_rank <<std::endl;
                std::fill(&input[local_length - h], &input[local_length], 0);
            }
            else { //Mixed case
                //std::cout << "Entering mixed case on rank " << world_rank <<std::endl;

                if (world_rank != world_size - 1){
                    MPI_Status status;
                    int position;
                    int number_amount;
                    //std::cout << "Awaiting recieve " <<world_rank<< std::endl;
                    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    position = status.MPI_TAG;
                    MPI_Get_count(&status, MPI_INT, &number_amount);


                    MPI_Recv(&input[position-offsets_start[world_rank]],number_amount, MPI_INT, MPI_ANY_SOURCE, position, MPI_COMM_WORLD, &status);//Recieve start
                    //std::cout << "Recieved " <<world_rank<< std::endl;

                    //Fill the rest with zeroes
                    std::fill(&input[number_amount], &input[local_length], 0);
                }
                else
                {
                    std::fill(&input[local_length - h], &input[local_length], 0);
                }


            }
            // std::cout << "Result on rank: " << world_rank << std::endl;
            // print_int_array(input, local_length);

    }

bool all_singleton (std::vector<tuple_ISA> &input, MPI_Comm comm, const int world_rank,
             const int world_size, int local_length){

    MPI_Datatype MPI_TUPLE_ISA;
    int lengths_ISA[2] = {1, 1};
    MPI_Aint displacements_ISA[2] = {0, sizeof(int)};
    MPI_Datatype types_ISA[2] = {MPI_INT, MPI_INT};
    MPI_Type_create_struct(2, lengths_ISA, displacements_ISA, types_ISA, &MPI_TUPLE_ISA);
    MPI_Type_commit(&MPI_TUPLE_ISA);

    if (world_rank != 0)
    {
        MPI_Send(&input[0], 1, MPI_TUPLE_ISA, world_rank - 1, 0, comm);
    }

    tuple_ISA after_end[1];
    if (world_rank != world_size - 1)
    {
        MPI_Recv(after_end, 1, MPI_TUPLE_ISA, world_rank + 1, 0, comm, MPI_STATUS_IGNORE);
    }
    else
    {
        after_end[0].SA = -1;
    }


    bool singleton = true;

    for (int i = 0; i < local_length - 1; i++){
        if (input[i].B == input[i + 1].B)
            singleton = false;
    }
    if (input[local_length - 1].B == after_end[0].B)
        singleton = false;

    //std::cout << singleton <<std::endl;
    return singleton;
}
