# BLAST_threading version 1.3
# Program for the automation of BLAST (Basic Local Alignment Search Tool) hosted by NCBI.
# Input is a file of txt format containing DNA sequences, output is data about de executed BLAST query from NCBI.
# Multiple queries are divided among different threads, the use of threads speeds up the progress significantly.
# Output data is being pushed to a mysql database.
# Last update: 10-06-2016

# Imports include modules of the Biopython package(version 5.5), mysqlconnector package(version..) and from the standard
# library: time and threading
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import mysql.connector
import time
import threading

# Log file is used to save data about the output of the program and any occurring errors
starttime = time.time()
log_file = open("Blast_log.txt", "a")
log_file.write(str(starttime))

# Contents of txt file are read and converted to a string, DNA sequences from txt file are subsequently converted
# to fasta format. NCBI BLAST queries are required to be fasta format. The first set of characters in the line
# containing a DNA sequence.

def fastaconvert(file):
    bestand = file.read()
    lines = bestand.replace("\t", "\n")
    lines = lines.split()

    dna = ["A", "T", "C", "G"]
    n_seq = []
    seq = []
    Header = False
    for line in lines:
        if "@HWI" in line:
            if Header:
                seq = "".join(seq)
                n_seq.append(seq)
                seq = []
            header = "".join(line.split(":")[5:])
            seq.append(">" + header + "\n")
            Header = True
        elif all(x in line for x in dna):
            seq.append(line)
    return n_seq

# Class with method to execute automated BLAST
# Objects consist of a name which is a number to annotate the individual threads
# The imported Biopython classes from the Bio.Blast module are used here: NCBIWWW and NCBIXML
# NCBIWWW is used to send the queries to NCBI, the parameters from qblast consist of: the type of BLAST, the sequence
# database with which the query sequences are being aligned with, the maximum number of results to be fetched,
# the e-value threshold, and the used score matrix

class BLAST_threads (threading.Thread):
    def __init__(self, name, seq_list):
        threading.Thread.__init__(self)
        self.name = name
        self.seq_list = seq_list
        
    def run(self):

        print("Starting thread " + self.name)
        seq_list = self.seq_list

        for seq in seq_list:
            notdone = True
            while notdone:
                try:
                    results_handle = NCBIWWW.qblast("blastx", "nr", seq, hitlist_size=5, expect=0.001, matrix_name="BLOSUM62")
                    blast_records = NCBIXML.parse(results_handle)
                    c = 1

                    # the number of hits are determined here
                    for blast_record in blast_records:
                        desc = blast_record.descriptions
                        numhits = len(desc)

                    # if the last number of the header of the sequence is 1: it is a forward read
                    # else (if the number of the header of the sequence is 2: it ia reverse read
                    header = seq.split()[0].lstrip(">")
                    sequentie = seq.split()[1]
                    print(header, " numhits: ", numhits)
                    if header.endswith("1"):
                        forward_reverse = "forward"
                    else:
                        forward_reverse = "reverse"
                    BLAST_hit = "1"
                    if numhits == 0:
                        BLAST_hit = "0"

                    # the mysql connector credentials are required to connect to the host containing the database
                    # the information about the dna sequence is pushed to the "sequentie_info" tabel of the db
                    conn = mysql.connector.connect(host="host", db="db", user="user",
                                                   passwd="passwd")
                    cursor = conn.cursor()
                    seq_info_query = ("INSERT INTO `pg2`.`sequentie_info` (`header`, `sequentie`, `forward_reverse`, `BLAST_hit`)\
                                                VALUES ('" + header + "', '" + sequentie + "', '" + forward_reverse + "', '" + BLAST_hit + "');")
                    cursor.execute(seq_info_query)
                    conn.commit()
                    cursor.close()
                    conn.close()

                    # boolean "notdone" is used to restart the while loop if a strange error occurres (e.g. error due to
                    # connection failure, other expected errors are being handled by excepts.

                    notdone = False
                except ValueError:
                    print("An error occurred during BLAST procedure")
                    log_file.write("An error occurred during BLAST procedure at: "+ str(header))
                    notdone = False
                except AttributeError:
                    print("An error occurred during BLAST procedure")
                    log_file.write("An error occurred during BLAST procedure at: " + str(header))
                    notdone = False
                except mysql.connector.ProgrammingError:
                    notdone = False
                    log_file.write("ProgrammingError occurred at: "+ str(header))
                except mysql.connector.DatabaseError:
                    notdone = False
                    log_file.write("DatabaseError occurred at: " + str(header))
                except mysql.connector.Error as err:
                    print("Oops, something went wrong: {}".format(err))
                    notdone = True

            # if the number of hits is not null, the information about the executed BLAST is fetched with the
            # hsp function, the data with tuple datatypes are subsequently converted to strings
            if numhits != 0:
                for alignment in blast_record.alignments:
                    notdone = True
                    while notdone:
                        try:
                            for hsp in alignment.hsps:
                                BLAST_hit = "1"
                                resultaat_id = str(header) + "(" + str(c) + ")"
                                accessie = str(alignment.accession)
                                hit_header = str(alignment.title)
                                hit_header = hit_header[:len(hit_header) + (100 - len(hit_header))]
                                score = str(hsp.score)
                                e_value = str(hsp.expect)
                                identities = str(hsp.identities)
                                positives = str(hsp.positives)
                                query_cover = str((float(identities) / float(len(hsp.query))) * 100)
                                gaps = str(hsp.gaps)
                                frame = str(hsp.frame)
                                BLAST_type = "blastx"

                                # the complete information about the BLAST is concatenated to a string for prints
                                # during the runtime of the program

                                alignmentnr = ("-----Alignment-" + str(c) + "----THREAD-" + self.name)
                                Blastresults = "\n" + alignmentnr + "\n" + "Hithead: " + hit_header + "header " + header + \
                                               "\n" \
                                               + "F/R: " + forward_reverse + " Accession: " + accessie + "\n" + "Score: " + str(
                                    score) \
                                               + " Identities: " + str(identities) + "\n" + "Positives: " + str(positives) \
                                               + " Gaps: " + str(gaps) + " Frame: " + str(frame) + "\n" + "Length: " + str(
                                    alignment.length) \
                                               + " E-value: " + str(e_value) + " BLAST_hit: " + str(BLAST_hit) \
                                               + "\n" + "resultaat_id: " + resultaat_id + " query cover: " + str(
                                    query_cover)

                                print(Blastresults)
                                log_file.write(Blastresults)

                                # connection to mysql db is reopened to push data about the executed BLAST
                                # when done, the connection is closed again
                                conn = mysql.connector.connect(host="host", db="db", user="user",
                                                               passwd="passwd")

                                cursor = conn.cursor()


                                BLAST_query = ("INSERT INTO `pg2`.`BLAST_resultaten` (`resultaat_id`, `accessie`, `hit_header`, `score`, `query_cover`,\
                                                          `e_value`, `identities`, `positives`, `gaps`, `frame`, `BLAST_type`, `sequentie_header`)\
                                                          VALUES ('" + resultaat_id + "', '" + accessie + "', '" + hit_header + "', '" + score + "', '" + query_cover + "', \
                                                            '" + e_value + "','" + identities + "', '" + positives + "','" + gaps + "', '" + frame + "', \
                                                            '" + BLAST_type + "', '" + header + "');")
                                cursor.execute(BLAST_query)
                                conn.commit()
                                cursor.close()
                                conn.close()
                                c += 1

                            # boolean "notdone" is used to restart the while loop if a strange error occurres (e.g. error due to
                            # connection failure, other expected errors are being handled by excepts.
                            notdone = False
                        except ValueError:
                            print("An error occurred during XML parsing")
                            log_file.write("An error occurred during XML parsing at: "+ str(resultaat_id))
                            notdone = False
                        except AttributeError:
                            print("An error occurred during XML parsing")
                            log_file.write("An error occurred during XML parsing at: "+ str(resultaat_id))
                            notdone = False
                        except mysql.connector.ProgrammingError:
                            notdone = False
                            log_file.write("ProgrammingError occurred at: "+ str(resultaat_id))
                        except mysql.connector.DatabaseError:
                            notdone = False
                            log_file.write("DatabaseError occurred at: "+ str(resultaat_id))
                        except mysql.connector.Error as err:
                            print("Oops, something went wrong: {}".format(err))
                            notdone = True

def main():


    # input file containing dna sequences is opened here in read mode, to prevent unwanted changes to the file.
    input_file = open("werkblad10.txt", "r")

    # get sequences in fastaformat from tab delineated txt file
    n_seq = fastaconvert(input_file)
    
    input_file.close()

    # four lists with unique data from n_seq, appended with steps of 10
    # e.g. list1 has index 0, 10, 20 etc, list2 has index 1, 11, 21 etc
    start = 0
    fasta_list1 = [n_seq[start] for start in range(start, int((len(n_seq)-(start)/10)), 10)]
    start = 1
    fasta_list2 = [n_seq[start] for start in range(start, int((len(n_seq)-(start)/10)), 10)]
    start = 2
    fasta_list3 = [n_seq[start] for start in range(start, int((len(n_seq)-(start)/10)), 10)]
    start = 3
    fasta_list4 = [n_seq[start] for start in range(start, int((len(n_seq)-(start)/10)), 10)]
    start = 4
    fasta_list5 = [n_seq[start] for start in range(start, int((len(n_seq)-(start)/10)), 10)]
    start = 5
    fasta_list6 = [n_seq[start] for start in range(start, int((len(n_seq)-(start)/10)), 10)]
    start = 6
    fasta_list7 = [n_seq[start] for start in range(start, int((len(n_seq)-(start)/10)), 10)]
    start = 7
    fasta_list8 = [n_seq[start] for start in range(start, int((len(n_seq)-(start)/10)), 10)]
    start = 8
    fasta_list9 = [n_seq[start] for start in range(start, int((len(n_seq)-(start)/10)), 10)]
    start = 9
    fasta_list10 = [n_seq[start] for start in range(start, int((len(n_seq)-(start)/10)), 10)]
    
    
    # Create new threads by creating objects in the BLAST_threads
    thread1 = BLAST_threads("1", fasta_list1)
    thread2 = BLAST_threads("2", fasta_list2)
    thread3 = BLAST_threads("3", fasta_list3)
    thread4 = BLAST_threads("4", fasta_list4)
    thread5 = BLAST_threads("5", fasta_list5)
    thread6 = BLAST_threads("6", fasta_list6)
    thread7 = BLAST_threads("7", fasta_list7)
    thread8 = BLAST_threads("8", fasta_list8)
    thread9 = BLAST_threads("9", fasta_list9)
    thread10 = BLAST_threads("10", fasta_list10)

    # Making a list with the above threads
    threads = []
    threads.extend([thread1, thread2, thread3, thread4, thread5, thread6, thread7, thread8, thread9, thread10])

    # Start new Threads
    for thread in threads:
        thread.start()
        time.sleep(3)
    # Wait for all threads to finish 
    for thread in threads:
        thread.join()

    # Runtime is documented in log file
    print("Exiting Main Thread")
    endtime = "It took", time.time()-starttime, "seconds"
    print(endtime)
    log_file.write("\n"+ str(endtime))
    log_file.close()

main()


