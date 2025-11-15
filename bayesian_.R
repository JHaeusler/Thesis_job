ggplot(df_plot, aes(x = p, y = Density, group = Scenario)) +
  
  # Mapear el color al factor 'Group' (1 a 6)
  geom_line(aes(color = factor(Group), linetype = Type), linewidth = 1) +
  
  # Usar linetypes para distinguir Uniforme de Beta dentro del mismo color
  scale_linetype_manual(values = c("Beta" = "dashed", "Uniforme" = "solid"), name = "Tipo Distr.") +
  
  # Usar una paleta de colores para los 6 grupos
  scale_color_manual(values = palette, name = "Grupo Escenario") +
  
  # Limitar el eje Y para mejor visualización
  ylim(0, y_limit) +
  xlim(0, 1) +
  # Título y etiquetas
  labs(title = "Curvas de Densidad A Priori: Uniforme vs. Beta (Agrupadas por Escenario)",
       subtitle = paste0("Mismo color = Mismo rango (Uniforme línea discontinua; Beta línea continua)"),
       x = "Proporción Defectuosa (p)",
       y = paste0("Densidad de Probabilidad (f(p)) - Máx ", y_limit)) +
  
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "right",
        legend.title = element_text(face = "bold"))