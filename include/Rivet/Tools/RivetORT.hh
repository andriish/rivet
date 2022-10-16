// -*- C++ -*-
#ifndef RIVET_RivetORT_HH
#define RIVET_RivetORT_HH

#include <iostream>

#include "Rivet/Tools/RivetPaths.hh"
#include "onnxruntime_cxx_api.h"



namespace Rivet {
  using namespace std;
    //Simple object to automatically take care of basic ONNX networks
    //Assumes one input/output node (note a node is not a neuron - a node is a single 
    //tensor of arbitrary dimension size)
    class RivetORT{

      public:
      RivetORT(const string& filename, const string& runname = "RivetORT"){
        //Set some ORT variables that need to be kept in memory
        _env = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING,runname.c_str() );
        Ort::SessionOptions sessionopts; //todo - check this is allowed to go out of scope.
        _session = std::make_unique<Ort::Session> (*_env, filename.c_str(), sessionopts);

        //Get netowrk hyper-params and store them (input, output shape, etc.)
        getNetworkInfo();

      }

      void compute(std::vector<float> &inputs, std::vector<float>& outputs){
        //Create ORT inputs
        // create input tensor object from data values
        auto memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
        auto input_tensor = Ort::Value::CreateTensor<float>(memory_info, inputs.data(),
                                                              inputs.size(), _inputNodeDims.data(), 2);//TODO: Understand this magic number

        //Work around for stupid pointer stuff
        const char* temp_inputNodeName = _inputNodeName.c_str();
        const char* temp_outputNodeName = _outputNodeName.c_str();
        
        auto output_tensors =
          _session->Run(Ort::RunOptions{nullptr}, &temp_inputNodeName, &input_tensor, 
                              1, &temp_outputNodeName, 1);//"magic" 1's reflect the number of input/output nodes.

        float* floatarr = output_tensors.front().GetTensorMutableData<float>();

        outputs.clear();
        //TODO: Generalise for different shape output arrays.
        for (size_t i = 0; i < _outputNodeDims[0]; ++i){
          outputs.push_back(*floatarr);
          ++floatarr;
        }                                                       

      }

      void getNetworkInfo(){
        Ort::AllocatorWithDefaultOptions allocator;
        //Magic 0's are the fact we only have 1 input/output node
        auto input_name = _session->GetInputNameAllocated(0, allocator);
        _inputNodeName = input_name.get();
        auto in_type_info = _session->GetInputTypeInfo(0);
        auto in_tensor_info = in_type_info.GetTensorTypeAndShapeInfo();
        _inType = in_tensor_info.GetElementType();//TODO: Use this for SFINAE
        _inputNodeDims = in_tensor_info.GetShape();

        auto output_name = _session->GetOutputNameAllocated(0, allocator);
        _outputNodeName = output_name.get();
        auto out_type_info = _session->GetInputTypeInfo(0);
        auto out_tensor_info = in_type_info.GetTensorTypeAndShapeInfo();
        _outType = out_tensor_info.GetElementType();//TODO: Use this for SFINAE
        _outputNodeDims = out_tensor_info.GetShape();

        //Do some basic sanity checks:
        if (_session->GetInputCount() != 1 || _session->GetInputCount() != 1){
          throw("RivetORT class cannot deal with multiple input/output nodes");
        }
      }

      std::ostream& printNetworkInfo(std::ostream& os){
          os << "RivetORT Network Summary: \n";
          os << "Input name: " << _inputNodeName << " Output name: " << _outputNodeName;
          os << "\nInput dimensions: (";
          for (size_t i = 0; i < _inputNodeDims.size()-1; ++i){
            os << _inputNodeDims[i] << ", ";
          }
          os << _inputNodeDims[_inputNodeDims.size() - 1] << ")\n";
          os << "Input Type (in ONNX enum form): " << _inType << "\n";
          os << "\nOutput dimensions: (";
          for (size_t i = 0; i < _outputNodeDims.size()-1; ++i){
            os << _outputNodeDims[i] << ", ";
          }
          os << _outputNodeDims[_outputNodeDims.size() - 1] << ")\n";
          os << "Output Type (in ONNX enum form): " << _outType << "\n";
          return os;
      }



      private:
      //ORT things that need to be preserved
      std::unique_ptr<Ort::Env> _env;
      std::unique_ptr<Ort::Session> _session;

      //Useful Info about the object, may as well be cached in member vars, will be needed
      //multiple times.
      std::vector<int64_t> _inputNodeDims;//I don't like this int64_t stuff but ORT insisted.
      std::vector<int64_t> _outputNodeDims;
      ONNXTensorElementDataType _inType;
      ONNXTensorElementDataType _outType;

      string _inputNodeName; //Needs to be converted for use but I'm not messing around with pointers to pointers to chars etc. Its 2022, guys!
      string _outputNodeName;

      

        

  };


}


#endif