# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import sys
import shutil
import hashlib
import subprocess
import requests
import re
import traceback
import uuid
from datetime import datetime
from pprint import pprint, pformat

import numpy as np
from Bio import SeqIO

from biokbase.workspace.client import Workspace as workspaceService
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from SetAPI.SetAPIClient import SetAPI
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from KBaseReport.KBaseReportClient import KBaseReport

#END_HEADER


class MEGAHIT:
    '''
    Module Name:
    MEGAHIT

    Module Description:
    A KBase module to wrap the MEGAHIT package.
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    VERSION = "3.0.0"
    GIT_URL = "https://github.com/dcchivian/kb_megahit"
    GIT_COMMIT_HASH = "fba68b5090e01bec120218e374f1b46fed077b4e"
    
    #BEGIN_CLASS_HEADER
    MEGAHIT = '/kb/module/megahit/megahit'

    # target is a list for collecting log messages
    def log(self, target, message):
        # we should do something better here...
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    # run megahit on the fastq files and return contig file path
    def exec_megahit_single_library (self, params):

        print('RUNNING MEGAHIT')
        print('Input reads files:')
        #fwd = reads[input_ref]['files']['fwd']
        #rev = reads[input_ref]['files']['rev']
        fwd = params['input_fwd_path']
        rev = params['input_rev_path']
        pprint('forward: '+fwd)
        pprint('reverse: '+rev)


        ### STEP 4: run megahit
        # construct the command
        megahit_cmd = [self.MEGAHIT]

        # we only support PE reads, so add that
        megahit_cmd.append('-1')
        megahit_cmd.append(fwd)
        megahit_cmd.append('-2')
        megahit_cmd.append(rev)

        # if a preset is defined, use that:
        if 'megahit_parameter_preset' in params:
            if params['megahit_parameter_preset']:
                megahit_cmd.append('--presets')
                megahit_cmd.append(params['megahit_parameter_preset'])

        if 'min_count' in params:
            if params['min_count']:
                megahit_cmd.append('--min-count')
                megahit_cmd.append(str(params['min_count']))
        if 'k_min' in params:
            if params['k_min']:
                megahit_cmd.append('--k-min')
                megahit_cmd.append(str(params['k_min']))
        if 'k_max' in params:
            if params['k_max']:
                megahit_cmd.append('--k-max')
                megahit_cmd.append(str(params['k_max']))
        if 'k_step' in params:
            if params['k_step']:
                megahit_cmd.append('--k-step')
                megahit_cmd.append(str(params['k_step']))
        if 'k_list' in params:
            if params['k_list']:
                k_list = []
                for k_val in params['k_list']:
                    k_list.append(str(k_val))
                megahit_cmd.append('--k-list')
                megahit_cmd.append(','.join(k_list))
        if 'min_contig_len' in params:
            if params['min_contig_len']:
                megahit_cmd.append('--min-contig-len')
                megahit_cmd.append(str(params['min_contig_len']))

        # set the output location
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        megahit_cmd.append('-o')
        megahit_cmd.append(output_dir)

        # run megahit
        print('running megahit:')
        print('    '+' '.join(megahit_cmd))
        p = subprocess.Popen(megahit_cmd, cwd=self.scratch, shell=False)
        retcode = p.wait()

        print('Return code: ' + str(retcode))
        if p.returncode != 0:
            raise ValueError('Error running MEGAHIT, return code: ' +
                             str(retcode) + '\n')

        output_contigs = os.path.join(output_dir, 'final.contigs.fa')
        if self.mac_mode: # on macs, we cannot run megahit in the shared host scratch space, so we need to move the file there
            shutil.move(output_contigs, os.path.join(self.host_scratch, 'final.contigs.fa'))
            output_contigs = os.path.join(self.host_scratch, 'final.contigs.fa')

        # send back path to contigs fasta file
        return output_contigs

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.serviceWizardURL = config['service-wizard-url']
        self.callbackURL = os.environ['SDK_CALLBACK_URL']
        self.scratch = os.path.abspath(config['scratch'])

        pprint(config)

        # HACK!! for issue where megahit fails on mac running docker because of 
        # silent named pipe error in the volume mounted from mac
        self.mac_mode = False
        self.host_scratch = self.scratch
        if 'mac-test-mode' in config and config['mac-test-mode']=='1':
            print('WARNING! running in mac test mode')
            self.scratch = os.path.join('/kb','module','local_scratch')
            self.mac_mode = True
        # end hack
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)
        #END_CONSTRUCTOR
        pass
    

    def run_megahit(self, ctx, params):
        """
        :param params: instance of type "MegaHitParams" (run_megahit() Run
           MEGAHIT.  Most parameters here are just passed forward to MEGAHIT
           run_megahit() is responsible for accepting input params from
           Narrative, calling exec_megahit(), and generating report. It
           mediates communication with the Narrative workspace_name - the
           name of the workspace for input/output input_reads_name - the name
           of the PE read library or ReadsSet (SE library support in the
           future) output_contig_set_name - the base name of the output
           contigset or AssemblySet combined_assembly_flag - if input is a
           ReadsSet, indicate combined Assembly megahit_parameter_preset -
           override a group of parameters; possible values: meta           
           '--min-count 2 --k-list 21,41,61,81,99' (generic metagenomes,
           default) meta-sensitive  '--min-count 2 --k-list
           21,31,41,51,61,71,81,91,99' (more sensitive but slower) meta-large
           '--min-count 2 --k-list 27,37,47,57,67,77,87' (large & complex
           metagenomes, like soil) bulk            '--min-count 3 --k-list
           31,51,71,91,99 --no-mercy' (experimental, standard bulk sequencing
           with >= 30x depth) single-cell     '--min-count 3 --k-list
           21,33,55,77,99,121 --merge_level 20,0.96' (experimental, single
           cell data) min_count - minimum multiplicity for filtering
           (k_min+1)-mers, default 2 min_k - minimum kmer size (<= 127), must
           be odd number, default 21 max_k - maximum kmer size (<= 127), must
           be odd number, default 99 k_step - increment of kmer size of each
           iteration (<= 28), must be even number, default 10 k_list - list
           of kmer size (all must be odd, in the range 15-127, increment <=
           28); overrides '--k-min', '--k-max', and '--k-step'
           min_contig_length - minimum length of contigs to output, default
           200 @optional megahit_parameter_preset @optional min_count
           @optional k_min @optional k_max @optional k_step @optional k_list
           @optional min_contig_len) -> structure: parameter "workspace_name"
           of String, parameter "input_reads_name" of String, parameter
           "output_contigset_name" of String, parameter
           "combined_assembly_flag" of Long, parameter
           "megahit_parameter_preset" of String, parameter "min_count" of
           Long, parameter "k_min" of Long, parameter "k_max" of Long,
           parameter "k_step" of Long, parameter "k_list" of list of Long,
           parameter "min_contig_len" of Long
        :returns: instance of type "MegaHitOutput" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_megahit

        SERVICE_VER = 'dev'

        print('Running run_megahit with params=')
        pprint(params)


        ### STEP 1: basic parameter checks + parsing
        required_params = ['workspace_name',
                           'input_reads_name',
                           'output_contigset_name'
                          ]
        for required_param in required_params:
            if required_param not in params or params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")


        ### STEP 2: call exec_megahit() - input params are the same, so just pass through
        exec_megahit_output = self.exec_megahit (ctx, params)


        ### STEP 3: save the report
        reportObj = {
            'objects_created':[],
            'text_message':exec_megahit_output['report_text']
        }
        for obj_ref in exec_megahit_output['output_contigset_refs']:
            reportObj['objects_created'].append({'ref':obj_ref, 'description':'Assembled contigs'})

        report = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})


        ### STEP 4: contruct the output to send back
        output = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }

        #END run_megahit

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_megahit return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def exec_megahit(self, ctx, params):
        """
        :param params: instance of type "ExecMegaHitParams" (exec_megahit()
           Actual execution of MEGAHIT Accepts ReadsSet or a ReadsLibrary as
           Input Creates Assembly object(s) as output. Will eventually also
           create AssemblySet object if input is a ReadsSet and not running a
           combined assembly Other vars same as run_megahit()) -> structure:
           parameter "workspace_name" of String, parameter "input_reads_name"
           of String, parameter "output_contigset_name" of String, parameter
           "combined_assembly_flag" of Long, parameter
           "megahit_parameter_preset" of String, parameter "min_count" of
           Long, parameter "k_min" of Long, parameter "k_max" of Long,
           parameter "k_step" of Long, parameter "k_list" of list of Long,
           parameter "min_contig_len" of Long
        :returns: instance of type "ExecMegaHitOutput" -> structure:
           parameter "report_text" of String, parameter
           "output_contigset_ref" of list of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN exec_megahit

        SERVICE_VER = 'dev'

        print('Running exec_megahit with params=')
        pprint(params)

        ### STEP 0: init
        token = ctx['token']
        wsClient = workspaceService(self.workspaceURL, token=token)
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token


        ### STEP 1: basic parameter checks + parsing
        required_params = ['workspace_name',
                           'input_reads_name',
                           'output_contigset_name'
                          ]
        for required_param in required_params:
            if required_param not in params or params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")


        ### STEP 2: determine if input is a ReadsLibrary or ReadsSet
        input_reads_ref = params['workspace_name']+'/'+params['input_reads_name']
        try:
            [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

            input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':input_reads_ref}]})[0]
            input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_info[TYPE_I])  # remove trailing version

        except Exception as e:
            raise ValueError('Unable to get reads object from workspace: (' + input_reads_ref +')' + str(e))

        accepted_input_types = ["KBaseSets.ReadsSet", "KBaseFile.PairedEndLibrary" ]
        if input_reads_obj_type not in accepted_input_types:
            raise ValueError ("Input reads of type '"+input_reads_obj_type+"' not accepted.  Must be one of "+", ".join(accepted_input_types))

        if input_reads_obj_type == "KBaseSets.ReadsSet":
            required_param = 'combined_assembly_flag'
            if required_param not in params or params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")


        ### STEP 3: get the list of library references
        if input_reads_obj_type == "KBaseFile.PairedEndLibrary":
            readsSet_ref_list   = [input_reads_ref]
            readsSet_names_list = [params['input_reads_name']]
 
        elif input_reads_obj_type == "KBaseSets.ReadsSet":
            readsSet_ref_list   = []
            readsSet_names_list = []

            try:
                setAPI_Client = SetAPI (url=self.serviceWizardURL, token=ctx['token'])  # for dynamic service
                input_readsSet_obj = setAPI_Client.get_reads_set_v1 ({'ref':input_reads_ref,'include_item_info':1})

            except Exception as e:
                raise ValueError('SetAPI FAILURE: Unable to get read library set object from workspace: (' + str(input_reads_ref)+")\n" + str(e))
            for readsLibrary_obj in input_readsSet_obj['data']['items']:
                readsSet_ref_list.append(readsLibrary_obj['ref'])
                NAME_I = 1
                readsSet_names_list.append(readsLibrary_obj['info'][NAME_I])

        else:
            raise ValueError ("Input reads of type '"+input_reads_obj_type+"' not accepted.  Must be one of "+", ".join(accepted_input_types))


        ### STEP 4: If doing a combined assembly on a ReadsSet, download reads one at a time and combine
        if input_reads_obj_type == "KBaseSets.ReadsSet" and params['combined_assembly_flag'] != 0:
            timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
            input_dir = os.path.join(self.scratch,'input.'+str(timestamp))
            if self.mac_mode: # on macs, we cannot run megahit in the shared host scratch space, so we need to move the file there
                input_dir = os.path.join(self.host_scratch,'input.'+str(timestamp))

            # connect to ReadsUtils Client
            try:
                readsUtils_Client = ReadsUtils (url=self.callbackURL, token=ctx['token'])  # SDK local
            except:
                raise ValueError("Unable to get readsUtils_Client\n" + str(e))

            # start combined file
            write_buf_size = 65536
            combined_input_fwd_path = os.path.join (input_dir, 'input_reads_fwd.fastq')
            combined_input_rev_path = os.path.join (input_dir, 'input_reads_rev.fastq')
            combined_input_fwd_handle = open (combined_input_fwd_path, 'w', write_buf_size)
            combined_input_rev_handle = open (combined_input_rev_path, 'w', write_buf_size)

            # add libraries, one at a time
            for this_input_reads_ref in readsSet_ref_list:
                try:
                    readsLibrary = readsUtils_Client.download_reads ({'read_libraries': [this_input_reads_ref],
                                                                      'interleaved': 'false'
                                                                      })
                except Exception as e:
                    raise ValueError('Unable to get reads object from workspace: (' + this_input_reads_ref +")\n" + str(e))

                this_input_fwd_path = readsLibrary['files'][this_input_reads_ref]['files']['fwd']
                this_input_rev_path = readsLibrary['files'][this_input_reads_ref]['files']['rev']

                # append fwd
                this_input_path = this_input_fwd_path
                cat_file_handle = combined_input_fwd_handle
                with open (this_input_path, 'r', read_buf_size) as this_input_handle:
                    while True:
                        read_data = this_input_handle.read(buf_size)
                        if read_data:
                            cat_file_handle.write(read_data)
                        else:
                            break
                os.remove (this_input_path)

                # append rev
                this_input_path = this_input_rev_path
                cat_file_handle = combined_input_rev_handle
                with open (this_input_path, 'r', read_buf_size) as this_input_handle:
                    while True:
                        read_data = this_input_handle.read(read_buf_size)
                        if read_data:
                            cat_file_handle.write(read_data)
                        else:
                            break
                os.remove (this_input_path)

            combined_input_fwd_handle.close()
            combined_input_rev_handle.close()


        ### STEP 5: finally run MEGAHIT
        exec_megahit_single_library_params = params
        output_assemblyset_contigset_paths = []
        output_contigset_path = None


        # PairedEndLibrary
        if input_reads_obj_type == "KBaseFile.PairedEndLibrary":
            try:
                readsUtils_Client = ReadsUtils (url=self.callbackURL, token=ctx['token'])  # SDK local
                readsLibrary = readsUtils_Client.download_reads ({'read_libraries': input_reads_ref,
                                                                  'interleaved': 'false'
                                                                  })
            except Exception as e:
                raise ValueError('Unable to get reads object from workspace: (' + input_reads_ref +")\n" + str(e))

            input_fwd_path = readsLibrary['files'][input_reads_ref]['files']['fwd']
            input_rev_path = readsLibrary['files'][input_reads_ref]['files']['rev']
            exec_megahit_single_library_params['input_fwd_file_path'] = input_fwd_path
            exec_megahit_single_library_params['input_rev_file_path'] = input_rev_path

            # the key line
            output_contigset_path = self.exec_megahit_single_library (exec_megahit_single_library_params)
            output_assemblyset_contigset_paths.append (output_contigset_path)
            
            os.remove (input_fwd_path) # files can be really big
            os.remove (input_rev_path)

        # ReadsSet combined (already downloaded and combined fastqs)
        elif input_reads_obj_type == "KBaseSets.ReadsSet" and params['combined_assembly_flag'] != 0:

            input_fwd_path = combined_input_fwd_path
            input_rev_path = combined_input_rev_path
            exec_megahit_single_library_params['input_fwd_file_path'] = input_fwd_path
            exec_megahit_single_library_params['input_rev_file_path'] = input_rev_path

            # the key line
            output_contigset_path = self.exec_megahit_single_library (exec_megahit_single_library_params)
            output_assemblyset_contigset_paths.append (output_contigset_path)

            os.remove (input_fwd_path) # files can be really big
            os.remove (input_rev_path)

        # ReadsSet uncombined (still have to download)
        elif input_reads_obj_type == "KBaseSets.ReadsSet" and params['combined_assembly_flag'] == 0:
            # connect to ReadsUtils Client
            try:
                readsUtils_Client = ReadsUtils (url=self.callbackURL, token=ctx['token'])  # SDK local
            except:
                raise ValueError("Unable to get readsUtils_Client\n" + str(e))

            # get libraries, one at a time, and run MEGAHIT
            output_assemblyset_contigset_paths = []
            for this_input_reads_ref in readsSet_ref_list:
                try:
                    readsLibrary = readsUtils_Client.download_reads ({'read_libraries': [this_input_reads_ref],
                                                                      'interleaved': 'false'
                                                                      })
                except Exception as e:
                    raise ValueError('Unable to get reads object from workspace: (' + this_input_reads_ref +")\n" + str(e))

                this_input_fwd_path = readsLibrary['files'][this_input_reads_ref]['files']['fwd']
                this_input_rev_path = readsLibrary['files'][this_input_reads_ref]['files']['rev']
                exec_megahit_single_library_params['input_fwd_file_path'] = this_input_fwd_path
                exec_megahit_single_library_params['input_rev_file_path'] = this_input_rev_path

                # the key line
                this_output_contigset_path = self.exec_megahit_single_library (exec_megahit_single_library_params)
                output_assemblyset_contigset_paths.append (this_output_contigset_path)
                
                os.remove (this_input_fwd_path) # files can be really big
                os.remove (this_input_rev_path)

        # just in case we've confused ourselves
        else:  
            raise ValueError ("error in logic")


        ### STEP 6: save the resulting assembly
        assemblyUtil = AssemblyUtil(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
        output_contigset_refs = []
        output_contigset_names = []
        for i,this_output_contigset_path in enumerate(output_assemblyset_contigset_paths):
            if len(output_contigset_paths) == 1:
                assembly_name = params['output_contigset_name']
            else:
                assembly_name = readsSet_names_list[i]+"_MEGAHIT_Contigs"

            this_output_data_ref = assemblyUtil.save_assembly_from_fasta({ 
                                                    'file':{'path':this_output_contigset_path},
                                                    'workspace_name':params['workspace_name'],
                                                    'assembly_name':assembly_name
                                                    })

            output_contigset_refs.append (this_output_data_ref)
            output_contigset_names.append (assembly_name)


        ### STEP 7: generate the report text

        # compute a simple contig length distribution for the report
        report = ''
        for i,this_output_contigset_path in enumerate(output_assemblyset_contigset_paths):

            report += "MEGAHIT run for Read Library: "+readSet_names_list[i]+"\n"
            report += "-------------------------------------------------------------\n"
            report += "\n"
            lengths = []
            for seq_record in SeqIO.parse(output_contigset, 'fasta'):
                lengths.append(len(seq_record.seq))

                report += 'ContigSet saved to: '+params['workspace_name']+'/'+output_contigset_names[i]+'\n'
                report += 'Assembled into '+str(len(lengths)) + ' contigs.\n'
                report += 'Avg Length: '+str(sum(lengths)/float(len(lengths))) + ' bp.\n'

                bins = 10
                counts, edges = np.histogram(lengths, bins)
                report += 'Contig Length Distribution (# of contigs -- min to max basepairs):\n'
                for c in range(bins):
                    report += '   '+str(counts[c]) + '\t--\t' + str(edges[c]) + ' to ' + str(edges[c+1]) + ' bp\n'


        ### STEP 8: contruct the output to send back
        output = { 'report_text': report,
                   'output_contigset_refs': output_contigset_refs
                 }

        #END exec_megahit

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method exec_megahit return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
