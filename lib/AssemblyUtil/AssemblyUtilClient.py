# -*- coding: utf-8 -*-
############################################################
#
# Autogenerated by the KBase type compiler -
# any changes made here will be overwritten
#
############################################################

from __future__ import print_function
# the following is a hack to get the baseclient to import whether we're in a
# package or not. This makes pep8 unhappy hence the annotations.
try:
    # baseclient and this client are in a package
    from .baseclient import BaseClient as _BaseClient  # @UnusedImport
except:
    # no they aren't
    from baseclient import BaseClient as _BaseClient  # @Reimport
import time


class AssemblyUtil(object):

    def __init__(
            self, url=None, timeout=30 * 60, user_id=None,
            password=None, token=None, ignore_authrc=False,
            trust_all_ssl_certificates=False,
            auth_svc='https://kbase.us/services/authorization/Sessions/Login',
            service_ver='release',
            async_job_check_time_ms=100, async_job_check_time_scale_percent=150, 
            async_job_check_max_time_ms=300000):
        if url is None:
            raise ValueError('A url is required')
        self._service_ver = service_ver
        self._client = _BaseClient(
            url, timeout=timeout, user_id=user_id, password=password,
            token=token, ignore_authrc=ignore_authrc,
            trust_all_ssl_certificates=trust_all_ssl_certificates,
            auth_svc=auth_svc,
            async_job_check_time_ms=async_job_check_time_ms,
            async_job_check_time_scale_percent=async_job_check_time_scale_percent,
            async_job_check_max_time_ms=async_job_check_max_time_ms)

    def _check_job(self, job_id):
        return self._client._check_job('AssemblyUtil', job_id)

    def _get_assembly_as_fasta_submit(self, params, context=None):
        return self._client._submit_job(
             'AssemblyUtil.get_assembly_as_fasta', [params],
             self._service_ver, context)

    def get_assembly_as_fasta(self, params, context=None):
        """
        Given a reference to an Assembly (or legacy ContigSet data object), along with a set of options,
        construct a local Fasta file with the sequence data.  If filename is set, attempt to save to the
        specified filename.  Otherwise, a random name will be generated.
        :param params: instance of type "GetAssemblyParams" (@optional
           filename) -> structure: parameter "ref" of String, parameter
           "filename" of String
        :returns: instance of type "FastaAssemblyFile" -> structure:
           parameter "path" of String
        """
        job_id = self._get_assembly_as_fasta_submit(params, context)
        async_job_check_time = self._client.async_job_check_time
        while True:
            time.sleep(async_job_check_time)
            async_job_check_time = (async_job_check_time *
                self._client.async_job_check_time_scale_percent / 100.0)
            if async_job_check_time > self._client.async_job_check_max_time:
                async_job_check_time = self._client.async_job_check_max_time
            job_state = self._check_job(job_id)
            if job_state['finished']:
                return job_state['result'][0]

    def _export_assembly_as_fasta_submit(self, params, context=None):
        return self._client._submit_job(
             'AssemblyUtil.export_assembly_as_fasta', [params],
             self._service_ver, context)

    def export_assembly_as_fasta(self, params, context=None):
        """
        A method designed especially for download, this calls 'get_assembly_as_fasta' to do
        the work, but then packages the output with WS provenance and object info into
        a zip file and saves to shock.
        :param params: instance of type "ExportParams" -> structure:
           parameter "input_ref" of String
        :returns: instance of type "ExportOutput" -> structure: parameter
           "shock_id" of String
        """
        job_id = self._export_assembly_as_fasta_submit(params, context)
        async_job_check_time = self._client.async_job_check_time
        while True:
            time.sleep(async_job_check_time)
            async_job_check_time = (async_job_check_time *
                self._client.async_job_check_time_scale_percent / 100.0)
            if async_job_check_time > self._client.async_job_check_max_time:
                async_job_check_time = self._client.async_job_check_max_time
            job_state = self._check_job(job_id)
            if job_state['finished']:
                return job_state['result'][0]

    def _save_assembly_from_fasta_submit(self, params, context=None):
        return self._client._submit_job(
             'AssemblyUtil.save_assembly_from_fasta', [params],
             self._service_ver, context)

    def save_assembly_from_fasta(self, params, context=None):
        """
        WARNING: has the side effect of moving the file to a temporary staging directory, because the upload
        script for assemblies currently requires a working directory, not a specific file.  It will attempt
        to upload everything in that directory.  This will move the file back to the original location, but
        if you are trying to keep an open file handle or are trying to do things concurrently to that file,
        this will break.  So this method is certainly NOT thread safe on the input file.
        :param params: instance of type "SaveAssemblyParams" (Options
           supported: file / shock_id / ftp_url - mutualy exclusive
           parameters pointing to file content workspace_name - target
           workspace assembly_name - target object name Uploader options not
           yet supported taxon_reference: The ws reference the assembly
           points to.  (Optional) source: The source of the data (Ex: Refseq)
           date_string: Date (or date range) associated with data. (Optional)
           contig_information_dict: A mapping that has is_circular and
           description information (Optional)) -> structure: parameter "file"
           of type "FastaAssemblyFile" -> structure: parameter "path" of
           String, parameter "shock_id" of type "ShockNodeId", parameter
           "ftp_url" of String, parameter "workspace_name" of String,
           parameter "assembly_name" of String
        :returns: instance of String
        """
        job_id = self._save_assembly_from_fasta_submit(params, context)
        async_job_check_time = self._client.async_job_check_time
        while True:
            time.sleep(async_job_check_time)
            async_job_check_time = (async_job_check_time *
                self._client.async_job_check_time_scale_percent / 100.0)
            if async_job_check_time > self._client.async_job_check_max_time:
                async_job_check_time = self._client.async_job_check_max_time
            job_state = self._check_job(job_id)
            if job_state['finished']:
                return job_state['result'][0]

    def status(self, context=None):
        job_id = self._client._submit_job('AssemblyUtil.status', 
            [], self._service_ver, context)
        async_job_check_time = self._client.async_job_check_time
        while True:
            time.sleep(async_job_check_time)
            async_job_check_time = (async_job_check_time *
                self._client.async_job_check_time_scale_percent / 100.0)
            if async_job_check_time > self._client.async_job_check_max_time:
                async_job_check_time = self._client.async_job_check_max_time
            job_state = self._check_job(job_id)
            if job_state['finished']:
                return job_state['result'][0]
