from __future__ import division
# misc routines for any command_line method
import sys,os
import libtbx.phil
from libtbx.utils import Sorry
import libtbx.phil.command_line
from phenix.utilities.arg_display_methods import arg_display_methods
from phenix.autosol.init_utils import init_utils

class headers(arg_display_methods,init_utils):
  def get_help(self,command_name,master_params,summary):
    print summary
    master_params.format(python_object=
          master_params.fetch(sources=[]).extract()).show()

  def raise_missing(self,what):
      raise Sorry("""\
          Missing file name for %(what)s :
          Please add %(what)s=file_name
          to the command line to specify %(what)s .""" % vars())

  def get_keyword_table(self,args_read,out=sys.stdout):
    if not args_read: args_read=[]
    args=[]
    self.keyword_table=None
    for arg in args_read:
      if arg=='keyword_table': # special case
        self.keyword_table=True
        print >>out,"Keyword Table will be written"
      else:
        args.append(arg)
    return args

  def get_params(self,command_name,master_params,args,out=sys.stdout,
      formatted_params=None):
    """
    Returns:
     done,master_params,new_params,changed_params,help
    """

    if not isinstance(master_params,libtbx.phil.scope):
      master_params=libtbx.phil.parse(master_params)
    argument_interpreter = libtbx.phil.command_line.argument_interpreter(
      master_phil=master_params, home_scope=command_name)

    if "--show_defaults" in args:  # just show defaults and return
      master_params.format(python_object=
          master_params.fetch(sources=[]).extract()).show()
      return True,None,None,None,False

    help=False
    if args is not None  and "--help" in args:
       help=True
       new_args=[]
       for arg in args:
         if arg!='--help': new_args.append(arg)
       if new_args:
         self.get_keyword_help(command_name,master_params,args[1:],out=out)
         return False,None,None,None,False
       args=new_args


    # read eff file if any, otherwise phil_objecs=[]
    args,phil_objects=\
     self.get_eff_file_from_args(args=args,master_params=master_params,out=out)

    # read formatted parameters, if any
    if formatted_params:
      edited_formatted_params=self.replace_overall_scope(
          formatted_params,command_name)
      phil_formatted_params=libtbx.phil.parse(edited_formatted_params)
      phil_objects.append(phil_formatted_params)

    for arg in args:
      if os.path.isfile(arg): continue
      try:
        arg=self.add_quote_to_parameter_def(arg)
        command_line_params = argument_interpreter.process(
           arg=arg)
        phil_objects.append(command_line_params)
      except KeyboardInterrupt: raise
      except Exception:
        if not arg.find("=")>-1:
          raise Sorry("\nSorry, cannot find the file or keyword: %s\n" % arg)
        else:
          raise Sorry("\nSorry, cannot interpret the keyword: %s\n" % arg)

    params = master_params.fetch(sources=phil_objects).extract()
    changed_params = master_params.fetch_diff(sources=phil_objects).extract()
    return False,master_params,params,changed_params,help

  def replace_overall_scope(self,text,name):
    # replace text between start and first { with name
    t=text.split("{")
    t[0]=name
    return "{".join(t)
  def get_eff_file_from_args(self,args=[],master_params=None,out=sys.stdout):
    # find an eff file in args
    input_eff_file=None
    edited_args=[]
    for arg in args:
      if input_eff_file is None and self.is_eff_file(arg): # 2013-01-10 allow = in file name
         input_eff_file=arg
      else:
         edited_args.append(arg)

    eff_file_params=None
    if input_eff_file is not None:
      print >>out,"Parameters taken from:",input_eff_file
      eff_file_params=libtbx.phil.parse(open(input_eff_file).read())
      if eff_file_params is not None and \
        not isinstance(eff_file_params,libtbx.phil.scope):
          eff_file_params=master_params.format(python_object = eff_file_params)
    if eff_file_params:
      eff_file_phil_objects=[eff_file_params]
    else:
      eff_file_phil_objects=[]

    return edited_args,eff_file_phil_objects


  def is_eff_file(self,file): #return True if this is a valid file for PHIL
    if not file or not os.path.isfile(file):
      return False
    import libtbx.phil.command_line
    import iotbx
    import iotbx.phil
    text=open(file).read()
    if not text.find("}")>-1 or not text.find("{")>-1:
         return False
    else:
      try:
         master_params=iotbx.phil.parse(text)
         return True
      except KeyboardInterrupt:  self.stop_after_KeyboardInterrupt()
      except Exception:
        return False

  def remove_blank_args(self,args,remove_ignore_blanks=True,out=sys.stdout):
    new_args=[]
    for arg in args:
      if arg.find("=")==-1:
        if remove_ignore_blanks and arg=='ignore_blanks': continue
        new_args.append(arg)
      else:
        spl=arg.split("=")
        if len(spl)>1 and spl[1]:
          if remove_ignore_blanks and spl[0]=='ignore_blanks': continue
          new_args.append(arg)
        else:
          if hasattr(self,'verbose') and self.verbose:
            print >>out,"Skipping command '"+str(arg)+"'"
    return new_args

  def create_local_dir(self,prefix='TEMP',start_number=0):
      import os
      # make a directory that does not exist TEMPxx
      for i in xrange(start_number,100000):
        temp_dir=prefix+str(i)
        if not os.path.exists(temp_dir):
          try:
            os.mkdir(temp_dir)
            return temp_dir
          except Exception: continue

      line="Failed to make TEMP directory..."
      raise AssertionError,line
