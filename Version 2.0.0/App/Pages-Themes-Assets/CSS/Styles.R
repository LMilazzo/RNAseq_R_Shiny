

#Call To in ui and add results to ui taglist
getStyles <- function(){
#Data Tables----

DataTableStyles <- "             
     table.dataTable {
        background-color: #1C1E22; /* Match medium gray background */
        border: 1px solid #677780; /* Subtle border */
        color: white; /* Light text */
        border-collapse: separate;
        border-spacing: 0;
     }
      
     table.dataTable thead th {
      background-color: #273238; /* Slightly darker for the header */
      color: #ffffff; /* White text for the header */
      text-align: center;
      padding: 8px;
      font-size: 16px ;
     }
     
     table.dataTable tbody td {
      background-color: #1C1E22;
      font-size: 14px;
      padding: 9px ;
     }
    
    table.dataTable th, table.dataTable td {
      border: 1px solid #677780; /* Subtle inner cell borders */
      padding: 6px; /* Comfortable padding */
    }

    .dataTables_wrapper .dataTables_filter input {
      border: 1px solid white;
      background-color: #4a4a4a;
      color: white;
      padding: 5px 10px;
      border-radius: 4px;
      width: 200px;
    }

    .dataTables_wrapper .dataTables_length select {
      border: 1px solid white;
      background-color: #4a4a4a;
      color: white;
      padding: 5px;
      border-radius: 4px;
    }
    
    .dataTables_wrapper .dataTables_length label {
      color: white;
    }

    .dataTables_wrapper .dataTables_paginate .paginate_button {
      background-color: #1C1E22;
      color: white !important;
      border: 1px solid white;
      font-size: 14px; /* Adjust text size */
      padding: 5px 8px;
      margin: 2px;
      border-radius: 4px;
    }

    .dataTables_wrapper .dataTables_paginate .paginate_button.current {
      background-color: #2aa4eb;
      color: black;
      font-weight: bold;
    }
    
    .dataTables_wrapper .dataTables_info {
      color: white;
    }

    .dataTables_wrapper .dataTables_paginate .previous 
    .dataTables_wrapper .dataTables_paginate .next {
      font-size: 14px;
      color: white !important; 
    }
      
    .dataTables_wrapper .dataTables_filter label {
      font-size: 14px;
      color: white;
    }
  
  "

custom_loader <- "
  .custom-loader {
        transform: rotateZ(45deg);
        perspective: 1000px;
        border-radius: 50%;
        width: 35px;
        height: 30px;
        color: #46a4fb;
        margin-left: 3px;
      }
        .custom-loader:before,
        .custom-loader:after {
          content: '';
          display: block;
          position: absolute;
          top: 0;
          left: 0;
          width: inherit;
          height: inherit;
          border-radius: 50%;
          transform: rotateX(70deg);
          animation: 1s spin linear infinite;
        }
        .custom-loader:after {
          color: #FF3D00;
          transform: rotateY(70deg);
          animation-delay: .4s;
        }

      @keyframes rotate {
        0% {
          transform: translate(-50%, -50%) rotateZ(0deg);
        }
        100% {
          transform: translate(-50%, -50%) rotateZ(360deg);
        }
      }

      @keyframes rotateccw {
        0% {
          transform: translate(-50%, -50%) rotate(0deg);
        }
        100% {
          transform: translate(-50%, -50%) rotate(-360deg);
        }
      }

      @keyframes spin {
        0%,
        100% {
          box-shadow: .2em 0px 0 0px currentcolor;
        }
        12% {
          box-shadow: .2em .2em 0 0 currentcolor;
        }
        25% {
          box-shadow: 0 .2em 0 0px currentcolor;
        }
        37% {
          box-shadow: -.2em .2em 0 0 currentcolor;
        }
        50% {
          box-shadow: -.2em 0 0 0 currentcolor;
        }
        62% {
          box-shadow: -.2em -.2em 0 0 currentcolor;
        }
        75% {
          box-shadow: 0px -.2em 0 0 currentcolor;
        }
        87% {
          box-shadow: .2em -.2em 0 0 currentcolor;
        }
      }
   
"


#Other----
Zoom <- "body {
  -moz-transform: scale(0.9, 0.9); /* Moz-browsers */
    zoom: 0.9; /* Other non-webkit browsers */
    zoom: 90%; /* Webkit browsers */
}"
#Tooltip----
tooltip <- "
.tip-mark {
      font-size: 25px;
      color: black;
      cursor: pointer;
      position: relative;
      display: inline-block;
      border: 1px solid black;
      border-radius: 100px;
      padding: 1px;
      margin: 0px;
      text-align: bottom;
    }

    .tool-tip {
      visibility: hidden;
      width: 200px;
      background-color: white;
      font-size: 12px;
      color: black;
      text-align: center;
      border-radius: 5px;
      padding: 0px;
      position: absolute;
      z-index: 999;
      bottom: 100%;
      left: 50%;
      margin-left: -100px;
      opacity: 0.6;
      transition: opacity 0.3s ease;
      border: 1px solid black;
    }

    .tip-mark:hover .tool-tip {
      visibility: visible;
      opacity: 1;
    }
"

collapsible <- "

.collapsible-container {
      width: 100%;
      max-width: 600px;
      margin: auto;
    }

    details {
      margin-bottom: 10px;
      background-color: #3e444c;
      border: 1px solid #272b30;
      border-radius: 15px;
    }
    
    summary {
      padding: 10px;
      font-size: 16px;
      font-weight: bold;
      background-color: #3e444c;
      border-radius: 15px;
      cursor: pointer;
      position: relative;
    }

    summary::after {
      content: '+';
      font-size: 20px;
      position: absolute;
      right: 10px;
      top: 30%;
      transform: translateY(-50%);
    }

    content: '-'; 
      details[open] summary::after {
    }

    summary:hover {
      background-color: #272b30;
    }

    .content {
      padding: 10px;
      background-color: #3e444c;
      border-radius: 15px;
    }

    /* Optional styles to add some interactivity */
    details[open] summary {
      background-color: #3e444c;
    }

"


#Return----
  return (
    css.Tags <-  tags$style(
                  HTML(DataTableStyles), 
                  HTML(Zoom), 
                  HTML(custom_loader),
                  HTML(tooltip),
                  HTML(collapsible)
                )
  )

}